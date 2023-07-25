"""
Cryogenic Thermometry Readout System
The implementation and testing modules for the cryogenic thermomtry readout system board to be implemented on the Simon's Observatory(SO) Large Aperture Telescope Receiver(LATR).

Developed By: N. Hosein, A. Manduca
University of Pennsylvania
"""

import os, sys, time, logging
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
from numpy.typing import ArrayLike
from labjack import ljm
from datetime import datetime
from scipy import signal, stats
from sys import getsizeof, stderr
from itertools import chain
from collections import deque
from reprlib import repr

matplotlib.use("tkagg")
logging.basicConfig(filename='coldtherm.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filemode='a')

class ColdThermometryReadout():
    def __init__(self, channels_to_read:str):
       """
        Readout class for the cold thermometry system.

        Parameters
        ----------
        channels_to_read : str
            Input string specifying the channels to be read. Input 'all' for all channels or the channels spearated by the
            commas: '40,41,42'

        Attributes
        ----------
        handle : LabJack Handle
            The LabJack handle object of the LabJack currently connected to the computer.

        setup_names : str
            List of strings to write(LabJack attribute).

        sample_rate : int
            The number of samples taken by the LabJack per stream per channel.

        scan_list : ArrayLike
            List of Modbus addresses to collect samples from, per scan(LabJack attribute).

        scan_amount : int
            The number of streams to be taken. Set to -1 for an infintie number of readings.

        channel_names : list
            List of channels being read.

        readoutDictionary : dict
            Dictionary storing the voltage(V), resistance(kOhms), temperatures(mK), and timestamps for each stream.
        
        readoutBufferSize : int
            Maximum size of `readoutDictionary` in bytes before it is cleared.

        stream_num : int
            The current stream number of the instance.

        upper_stage_config : dict
            Dictionary storing the upper stage bias configuration of the board.

        middle_stage_config : dict
            Dictionary storing the middle stage bias configuration of the board.
            
        lower_stage_config : dict
            Dictionary storing the lower stage bias configuration of the board.
            
        ResistorCalibrationDictionary : dict
            Dictionary containing the calibration curve data for each resistor/channel.
            
        biasMode : int
            The DC signal in volts to be outputted from the output port set by `BiasOutputPort`.
            0: LS; 1: MS; 2: US; -1: Override

            Set to an abitrary integer that is not one of the above to set the system in automatic bias switching mode.

        BiasSwitchAveragingInterval : int
            The interval over which the board averages to determine the current temperature stage.

        overrideBias:float
            The manual voltage bias to be supplied to the signal generator chip if the automatic switching functionailty is no longer wanted.

        BiasOutputPort:str, default = 'DAC0'
            Name of the DAC port from which the DC voltage to generate the signal is sent.
        
        ResistorCalibrationMappingPath : str
            String specifying the location of the csv file which contains the mapping to the calibration file of each resistor.
        
        ReferenceChannel : str
            The name of the drift reference channel. Eg. 'AIN84'

        ReferenceChannelResistance : float
            The value of the fixed resistances attached to the reference channel in kOhms.
        """
       
       #Configuable Settings
       self.sample_rate:int = 1600
       self.scan_amount:int = 10
       self.readoutBufferSize:int = 5000
       self.upper_stage_config:dict = {"UL":100,"LL":1.2,"Bias":5, "CalibratedBias":5}
       self.middle_stage_config:dict = {"UL":1.2,"LL":0.15,"Bias":0.2, "CalibratedBias":0.2}
       self.lower_stage_config:dict = {"UL":0.15,"LL":0,"Bias":0.02, "CalibratedBias":0.02}
       self.BiasSwitchAveragingInterval:int = 30
       self.overrideBias:float = 3
       self.BiasOutputPort:str = "DAC0"
       self.ResistorCalibrationMappingPath:str = "./resistor_calibration_mapping.csv"
       self.ReferenceChannel:str = "AIN84"
       self.ReferenceChannelResistance:float = 20
    
       #Non-Configurable Settings
       self.handle = None
       self.setup_names:str = ''
       self.setup_values:str = ''
       self.scan_list:list = None
       self.channel_names:list = []
       self.readoutDictionary:dict = {}
       self.stream_num:int = 1
       self.ResistorCalibrationDictionary = {}
       self.biasMode:int = 10 #Arbitrarily Set
       
       self.OpenConnection()
       self.SetupLabJack()
       self.SetChannelsToRead(channels_to_read)
       self.GenerateDictionaries()
       self.LoadResistorCalibration(self.ResistorCalibrationMappingPath)

       self.TestingModule = Testing(self)

    def ClearMemoryBuffer(self):
        """
        Checks whether to clear the memory buffer `readoutDictionary`. Buffer is cleared if the size of `readoutDictionary` in bytes is 
        greater than the value set in `readoutBufferSize`.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        print self.total_size(self.readoutDictionary)
        if self.total_size(self.readoutDictionary) >= self.readoutBufferSize:
            self.readoutDictionary.clear()
            self.GenerateDictionaries()

    def total_size(o, handlers={}, verbose=False):
        """ 
        Checks the total memory utilized by a given container in bytes. See source below for further usage information.

        Source: https://code.activestate.com/recipes/577504/
        """
        dict_handler = lambda d: chain.from_iterable(d.items())
        all_handlers = {tuple: iter,
                        list: iter,
                        deque: iter,
                        dict: dict_handler,
                        set: iter,
                        frozenset: iter,
                    }
        all_handlers.update(handlers)     # user handlers take precedence
        seen = set()                      # track which object id's have already been seen
        default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

        def sizeof(o):
            if id(o) in seen:       # do not double count the same object
                return 0
            seen.add(id(o))
            s = getsizeof(o, default_size)

            if verbose:
                print(s, type(o), repr(o), file=stderr)

            for typ, handler in all_handlers.items():
                if isinstance(o, typ):
                    s += sum(map(sizeof, handler(o)))
                    break
            return s

        return sizeof(o)

    def BiasSwitch(self, bias:float):
        """
        Changes the bias signal from the LabJack to the signal generator chip. This value in most cases will not correspond to the output of the signal generator chip.

        Parameters
        ----------
        bias : float
            The DC signal in volts to be outputted from the output port set by `BiasOutputPort`.

        Returns
        -------
        None
        """

        ljm.eWriteName(self.handle, self.BiasOutputPort, bias)

    def CheckBiasSwitch(self, ovrr:bool = False):
        """
        Checks whether the bias signal needs to be switched. By default the the bias signal is set automatically by the values set in the dictionaries `lower_stage_config`, `middle_stage_config`, and `upper_stage_config`.
        The interval over which the temperature values are averaged to look for stage changes are given by `BiasSwitchAveragingInterval`.

        Parameters
        ----------
        ovrr : bool, default = False
            Enable manual mode to override the automatic functionality and set the output bias to a fixed value given by `overrideBias`.

        Returns
        -------
        None
        """
      
        if ovrr and self.biasMode == -1:
            return
        elif ovrr:
            self.BiasSwitch(self.overrideBias)
            return

        lowestChannelTemp = 1000
        for channel in self.channel_names:
            
            channeltemp = np.average(self.readoutDictionary[channel]["Temp [mK]"][-self.BiasSwitchAveragingInterval:])/1000 #Temp in K
            if channeltemp < lowestChannelTemp:
                lowestChannelTemp = channeltemp

        if lowestChannelTemp >= self.upper_stage_config["LL"]:
            if self.biasMode != 0:
                self.BiasSwitch(self.upper_stage_config["CalibratedBias"])

        elif lowestChannelTemp < self.upper_stage_config["UL"] and lowestChannelTemp >= self.middle_stage_config["LL"]:
            if self.biasMode != 1:
                self.BiasSwitch(self.middle_stage_config["CalibratedBias"])

        elif lowestChannelTemp < self.lower_stage_config["UL"]:
            if self.biasMode != 2:
                self.BiasSwitch(self.lower_stage_config["CalibratedBias"])

    def LoadResistorCalibration(self, path:str="./resistor_calibration_mapping.csv"):
        """
        Loads the resistor calibration curve for each resistor from the calibration file mapping `resistor_calibraton_mapping.csv`.

        Parameters
        ----------
        path : str, default = "./resistor_calibration_mapping.csv"
            The path to the file `resistor_calibraton_mapping.csv` which contains the mappings to the calibration curves of each resistor attached to each channel.

        Returns
        -------
        None
        """

        try:  
            mapping_file = os.path.realpath(path)
            df = pd.read_csv(mapping_file, header=0,names=["Channel Name", "Path"])
            channel_names = df["Channel Name"].values
            calibration_curve_paths = df["Path"].values
        except Exception:
            error = sys.exc_info()[1]
            print("Could not load resistor-calibration mapping file. Check file path.")
            quit()
        
        mappingdictionary = {}

        #Loads mappings into a dictionary
        for i in range(len(channel_names)):
            mappingdictionary[channel_names[i]]  = {"Path":calibration_curve_paths[i]}
        
        for channel_name in self.channel_names: #only loads the channels which are in use.
            try:
                path = os.path.realpath(mappingdictionary[channel_name]["Path"])
             
                with open(path) as calibrationfile:
                    data = np.array([])
                    for line in calibrationfile:
                        data = np.append(data, np.array(line.split(), dtype=float)) #This generates an array of size (n*4)

                    data = np.resize(data, (int(len(data)/4),4))
                    temp = data[:,0]
                    res = data[:,1]
                    drdt = data[:,2]
                    sd = data[:,3]
                   
                    self.ResistorCalibrationDictionary[channel_name] = {"Temp [K]":temp[::-1], "log(R [Ohms])":np.log10(res[::-1]), "drdt": drdt[::-1], "sd": sd[::-1]} #reverses array for np.interp

            except Exception:
                error = sys.exc_info()[1]
                print("An error occured while loading the calibration for channel: " + channel_name)        
    
    def GenerateDictionaries(self):
        """
        Regenerates the `readoutDictionary` when the buffer is filled. Buffer size is set by `readoutBufferSize`.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.readoutDictionary.clear()

        for n in self.channel_names:
            self.readoutDictionary[n] = {"V [V]":[],"R [kohms]":[],"Temp [mK]":[],"Time":[]}

    def OpenConnection(self):
        """
        Opens a connection to a T7 LabJack and stores the LabJack handle as `handle`. Closes connections to other LabJacks before executing.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """
       
        try:
            ljm.closeAll()

            logging.info("Openning connection to LabJack.")

            self.handle = ljm.openS("T7", "ANY", "ANY")  

            # grab and print out important info 
            info = ljm.getHandleInfo(self.handle)
            print("\nOpened a LabJack with Device type: %i, Connection type: %i,\n"
                "Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" %
                (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))
            
            logging.info("Opened a LabJack with Device type: %i, Connection type: %i,"
                " Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" %
                (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))
       
        except ljm.LJMError:
            
            ljme = sys.exc_info()[1]
            if ljme.errorCode == 1314: #No Devices found
                logging.WARNING("LabJack Error 1314 - No Devices Found.")

                print("\n Unable to find labjack device. Retrying in 5s.")
                time.sleep(5)

                self.OpenConnection()
    
    def SetChannelsToRead(self, channels_to_read:str):
        """
        Sets the channels to be read by the LabJack.

        Parameters
        ----------
        channels_to_read : str
            Input string specifying the channels to be read. Input 'all' for all channels or the channels spearated by the
            commas: '40,41,42'

        Returns
        -------
        None
        """

        try:
            channel_names = []
            if channels_to_read == 'all':

                # analog input channels to be read out
                positive_channels = np.append(np.append(np.arange(48,56), np.arange(80,88)), np.arange(96,104))
                negative_channels = np.append(np.append(np.arange(56,64), np.arange(88,96)), np.arange(104,112))
                channel_numbers = [item for sublist in zip(positive_channels, negative_channels) for item in sublist]

                for c in channel_numbers: channel_names = channel_names+["AIN%d"%c]
                #print(channel_names)
            else:
                for i in channels_to_read.split(","):
                    i = "AIN"+i
                    channel_names.append(i)
                #print(channel_names)

            self.scan_list = ljm.namesToAddresses(len(channel_names), channel_names)[0]
            self.channel_names = channel_names

            logging.info("Readout Channels Set: " + channels_to_read)
        
        except Exception:
            print(sys.exc_info()[1])
            raise
    
    def SetupLabJack(self, setup_names:ArrayLike=None, setup_values:ArrayLike=None):
        """
        Sets the LabJack.

        Parameters
        ----------
        setup_names : ArrayLike
            Array of names to write. See LabJack documentation.

        setup_values : ArrayLike
            Array of values to write. See LabJack documentation.
        Returns
        -------
        None

        Notes
        -----
       
        set range, resolution, and settling time for channels factor
        note:
        Negative channel: single ended = 199, differential = 1
        Range: The instrumentation amplifier in the T7 provides 4 different gains:
        x1 (RANGE is ±10 volts), enter 10.0
        x10 (RANGE is ±1 volts), enter 1.0
        x100 (RANGE is ±0.1 volts), enter 0.1
        x1000 (RANGE is ±0.01 volts), enter 0.01
        Resolution index = Default (0)
        Settling, in microseconds = Auto (0) resource on settling times: https://old3.labjack.com/support/app-notes/SettlingTime
                
        """

        if setup_names == None:
            self.setup_names = ["AIN_ALL_NEGATIVE_CH", "AIN_ALL_RANGE","STREAM_RESOLUTION_INDEX","STREAM_SETTLING_US"]
        else:
            self.setup_names = setup_names

        if setup_values == None:
            self.setup_values = [199,10.0,0,0]
        else:
            self.setup_values = setup_values

        try:
        # assign the values of range, resolution, and settling to each channel
            ljm.eWriteNames(self.handle, len(self.setup_names), self.setup_names, self.setup_values)
        except ljm.LJMError:
            
            ljme = sys.exc_info()[1]
            if ljme.errorCode == 2605: #Stream is Active.
                print("\n Stream is active. Stopping stream.")
                ljm.eStreamStop(self.handle)
                print("Stream stopped. Restarting Method.")
                self.SetupLabJack(setup_names, setup_values)


    def ConvertResistance(self, voltage:ArrayLike):
        """
        Converts a given array of voltages into rsistances using a pre-determined calibration equation. See `CalibrateBoard()`of the Testing class for determining the 
        calibration equation.

        Parameters
        ----------
        voltage : ArrayLike
            ArrayLike object containing the voltage output of the board in volts.

        Returns
        -------
        res : ArrayLike
            ArrayLike object containing the resistances of the thermometer in kOhms.
        """

        res = 58.69299438396007*voltage + 0.12049023712211593

        return res

    def Stream(self):
        """
        Streams the readout channels set in `channel_names` for the number of streams specified in `scan_amount`.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        try:
            ljm.eStreamStart(self.handle, int(self.sample_rate), len(self.channel_names), self.scan_list, self.sample_rate)
        except:
            ljme = sys.exc_info()[1]
            if ljme.errorCode == 2605:
                print("\n Stream is active. Stopping stream.")
                ljm.eStreamStop(self.handle)
                print("Stream stopped. Restarting Method.")
                self.Stream()

        print("\nStarting %s read loops.\n" % (str(self.scan_amount)))

        total_scans = 0
        errorcount = 0

        starttime = datetime.now()
        while self.stream_num <= self.scan_amount:
            try: 
                
                streamtimeutc = time.time()
                v_measured = ljm.eStreamRead(self.handle)[0]

                total_scans += 1
                
                for k in range(len(self.channel_names)):
                    
                    v = np.array(v_measured[k::len(self.channel_names)])
                    res = np.array(self.ConvertResistance(v))
                    temp = np.array(self.ConvertTemperature((res), self.channel_names[k]))

                    self.readoutDictionary[self.channel_names[k]]["V [V]"].append(v)
                    self.readoutDictionary[self.channel_names[k]]["R [kohms]"].append(res)
                    self.readoutDictionary[self.channel_names[k]]["Temp [mK]"].append(temp)
                    self.readoutDictionary[self.channel_names[k]]["Time"].append(streamtimeutc)

                self.CheckBiasSwitch()

                self.TestingModule.TestingRoutine()

                self.ClearMemoryBuffer()

            except ljm.LJMError:
                errorcount += 1
                ljme = sys.exc_info()[1]
                
                if ljme.errorCode == 1301: #buffer full:
                    print("Buffer Full")
                    
                elif ljme.errorCode == 1263: #No bytes Received.
                    print("No bytes received. Restarting Stream.")

                elif ljme.errorCode == 1306: #Synch Timeout
                    print("Synch Timeout.")

                elif ljme.errorCode == 1224: #device not open
                    print("Device not open. Reopenning connection.")       

                else:
                    print("Unknown Error")

                time.sleep(1)
                self.OpenConnection()
                self.Stream()

                pass

            self.stream_num += 1

        ljm.eStreamStop(self.handle)
        endtime = datetime.now()
        processtime = endtime - starttime
        
        print("Streaming Complete. \nTotal Scans: {2}; Errors Encountered: {0}; Stream Time: {1}".format(str(errorcount), str(processtime), str(total_scans)))

    def ConvertTemperature(self, resistance:ArrayLike, channel_name:str, unit:str="mK"):
        """
        Converts a given resistance value to a temperature value based on the resistance-temperature calibration of the given resistor.

        Parameters
        ----------
        resistance : ArrayLike
            An array containing the the resistance values(in kOhms) to be converted into temperature values.

        channel_name : str
            The name of the channel from which the resistance values were obtained.

        unit : str, default = 'mK'
            The unit to which the final temperature values are saved. 
            'mK' for millikelvin; 'K' for kelvin.

        Returns
        -------
        temp : ArrayLike
            An array containing the converted temperature values for the given channel.
        """
        factor = 0
        if unit == "mK":
            factor = 1000
        if unit == "K":
            factor = 1
        else:
            pass

        temp = factor * np.interp(np.log10(1000 * resistance), self.ResistorCalibrationDictionary[channel_name]["log(R [Ohms])"], self.ResistorCalibrationDictionary[channel_name]["Temp [K]"])
        
        return temp

    def CloseConnection(self):
        """
        Closes the connection to the current LabJack handle

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        ljm.close(self.handle)
        print("Connection closed.")

    def StreamChannel(self, channel:str, numStreams:int):
        """
        Streams a specified channel for the given number of streams and returns an array of the streams. The sample rate is set by `sample_rate`.

        Parameters
        ----------
        channel : str
            Input string specifying the channel to be read. Eg: '40'

        numStreams : int
            The number of samples streams to save.

        Returns
        -------
        voltages : ArrayLike
            ArrayLike object the voltage output of the streams.
        """

        self.SetChannelsToRead(channel)
        voltages = []

        try:
            ljm.eStreamStart(self.handle, int(self.sample_rate), len(self.channel_names), self.scan_list, self.sample_rate)
        except ljm.LJMError:
            ljme = sys.exc_info()[1]
            if ljme.errorCode == 2605:
                ljm.eStreamStop(self.handle)
                self.StreamChannel(channel, numStreams)

        streamNum = 0
      
        while streamNum < numStreams:
            try: 
                v_measured = ljm.eStreamRead(self.handle)[0]
                for k in range(len(self.channel_names)):
                    voltages.append(v_measured[k::len(self.channel_names)])
            
            except ljm.LJMError:
                ljme = sys.exc_info()[1]
                pass

            streamNum += 1

        return voltages

class Testing():
    def __init__(self, readoutObject:ColdThermometryReadout = None):
        """
        Testing class for the cold thermometry system.

        Parameters
        ----------
        readoutObject : ColdThermometryReadout, default = None
            The readout object of the current readout instance. When loading and analysing data, set to `None`.

        Attributes
        ----------
        readoutObject : ColdThermometryReadout, default = None
            The readout object of the current readout instance. When loading and analysing data, set to `None`.

        avg_dictionary : Dictionary
            The dictionary which stores the averaged data of the current stream instance.

        bin_dictionary : Dictionary
            The dictionary which stores the binned data of the current stream instance.

        LoadedDataDictionaryAveraged : Dictionary
            The dictionary which stores the loaded averaged data.

        LoadedDataDictionaryBinned : Dictionary
            The dictionary which stores the loaded binned data.

        saveinterval : int, default = 30
            The interval over which data is saved and memory beffers are cleared. An interval of 30 means that data is saved after 30 streams which results in a save rate of 2 saves per minute.

        numBins : int, default = 10
            The number of bins that the data is split into. The default sample rate of the board is 1600Hz. `numBins` splits the data into 10 bins of length 1600/10 = 160.

        saveCounter : int, default = 0
            The counter which stores how many files have been saved for the current readout instacne.

        channels_to_print : ArrayLike
            An array containing the names of the channels to print.
            Eg. ["AIN50", "AIN51"]
        
        channels_to_plot : ArrayLike
            An array containing the names of the channels to show on the live readout.
            Eg. ["AIN50", "AIN51"]

        fig : Figure
            The figure of the live readout graph.

        ax : Axis
            The axis object of the live readout graph.
        """
        
        #Configurable Settings
        self.channels_to_print = ["AIN56"]
        self.channels_to_plot = ["AIN56"]
        self.saveInterval = 30
        self.numBins = 10
        self.savePath = "./Board Test/"

        #Non-Configurable Settings
        self.readout = readoutObject
        self.avg_dictionary = {}
        self.bin_dictionary = {}
        self.LoadedDataDictionaryAveraged = {}
        self.LoadedDataDictionaryBinned = {}
        self.saveCounter = 0
        self.fig = None
        self.ax = None

        self.GenerateDictionaries()
        self.setup_animation()

    def GenerateDictionaries(self):
        """
        Regenerates the `avg_dictionary` and `bin_dictionary` when the buffer is filled. Buffer size is set by `saveInterval`.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.readout == None:
            return

        self.avg_dictionary.clear()
        self.bin_dictionary.clear()

        for n in self.readout.channel_names:
            self.avg_dictionary[n] = {"V [V]":[],"R [kohms]":[],"Temp [mK]":[],"Time":[]}
            self.bin_dictionary[n] = {"V [V]":[],"R [kohms]":[],"Temp [mK]":[],"Time":[]}

    def CalibrateBoard(self, channel:str, streamAmount:int):
        """
        Interface to calibrate the board and obtain a relationship between voltage(V) and resistance(kOhms).

        Parameters
        ----------
        channel : str
            The name of the channel from which the voltage samples will be taken. Eg: "AIN56"

        streamAmount : int
            The number of streams for a given resistance. 

        Returns
        -------
        None
        """

        self.readout.SetChannelsToRead(str(channel))

        print("\n ### Calibration ###")
        print("\nFollow the following instructions to calibrate the board.")

        print("You will be repeatedly prompted to enter resistances. These are the input resistances to the set channel. The application will then stream the channel for the specified number of samples and aggregate all data to determine the calibration for the board.")

        print("\nNOTE: Enter all resistances in kOhms. Press 'c' at any time to stop and complete the complete the calibration.")

        inputResistances = []
        outputVoltageAverage = []
        count = 1

        res = input("\nData Point {0}: Input resistance to channel {1}: ".format(str(count), str(channel)))
     
        while res != 'c':
            inputResistances.append(float(res))
            print("  => Streaming data. Please wait.")
            
            voltage = self.readout.StreamChannel(channel, streamAmount)
            avgVoltage = np.average(voltage)
            outputVoltageAverage.append(avgVoltage)
        
            print("  => Average Voltage: {0}V".format(avgVoltage))
            
            count += 1
            res = input("\nData Point {0}: Input resistance to channel {1}: ".format(str(count), str(channel)))

        print("\n Data Summary: (Voltage, Resistance)")
        for i in range(len(inputResistances)):
            print("({0}, {1})".format(outputVoltageAverage[i], inputResistances[i]))

        slope, intercept, r, p, se = stats.linregress(outputVoltageAverage, inputResistances)

        print("\nCalibration Function: res = {0}*voltage + {1}".format(slope, intercept))
        print("r: {0}".format(r))

        print("\nEnter the above equation into the calibraton function of the readout class.")

        print("\n### Calibration Complete ###")

    def redraw(self):
        """
        Redraws the graph in the live-readout graph.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
    
    def setup_animation(self):
        """
        Sets the interactive property of the graph if a live readout is in progress, and assigns the figure and axes objects to the Testing object.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if self.readout == None:
            plt.ioff()
            self.fig = plt.figure() 
            self.ax = self.fig.subplots()
        else:
            plt.ion() 
            self.fig = plt.figure() 
            self.ax = self.fig.subplots() 

    def Plot(self, channels:ArrayLike, attr:str = "V [V]", graphtitle:str = "Cold Thermometry Test Readout"):
        """
        Defines which channels to plot for the live readout.

        Parameters
        ----------
        channels : Arraylike
            An array of the channels to be read. Each element of the array must be the channel name.
            Eg: ["AIN50", "AIN51", "AIN52"]
        
        attr : str, default = "V [V]"
            The attribute which will be plotted against time.
            For Voltage, attr = "V [V]";
            For Resistance, attr = "R [kohms]";
            For Temperature, attr = "Temp [mK]".

        graphtitle : str, default = "Cold Thermometry Test Readout"
            Sets the title of the graph in the live readout.

        Returns
        -------
        None
        """
            
        plt.cla()
        self.redraw()
            
        for channel in channels:
            self.ax.plot(self.avg_dictionary[channel][attr], label = str(channel))
        
        plt.title(graphtitle)
        plt.xlabel("t [s]")
        plt.ylabel(attr)
        plt.legend()
        plt.pause(0.5)

    def CalculateSplitTime(self):
        """
        Calculates the times for the binned arrays/dictionaries.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for channel in self.avg_dictionary:
            #Use onyly first channel as reference
            timediff = 1/self.numBins
            array = np.array([])

            for time in self.avg_dictionary[channel]["Time"]:
                for j in range(self.numBins):
                    array = np.append(array, time + j*timediff)

            return array

    def AppendAverageReadout(self):
        """
        Appends the average of the last stream of the current readout instance to the `avg_dictionary`.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for channel in self.readout.channel_names:
            
            self.avg_dictionary[channel]["V [V]"].append(np.average(self.readout.readoutDictionary[channel]['V [V]'][-1]))
            self.avg_dictionary[channel]["R [kohms]"].append(np.average(self.readout.readoutDictionary[channel]['R [kohms]'][-1]))
            self.avg_dictionary[channel]["Temp [mK]"].append(np.average(self.readout.readoutDictionary[channel]['Temp [mK]'][-1]))
            self.avg_dictionary[channel]["Time"].append(np.average(self.readout.readoutDictionary[channel]['Time'][-1])) 

    def AppendBinReadouts(self):
        """
        Appends the bin readouts to `bin_dictionary`. The number of bins is set in `numBins`.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for channel in self.readout.channel_names:

            voltage_arr = self.readout.readoutDictionary[channel]["V [V]"][-1]
            voltage_bins = np.reshape(voltage_arr, (self.numBins, int(len(voltage_arr)/self.numBins)))
            voltage_bins = np.average(voltage_bins, 1)

            res_arr = self.readout.readoutDictionary[channel]["R [kohms]"][-1]
            res_bins = np.reshape(res_arr, (self.numBins, int(len(voltage_arr)/self.numBins)))
            res_bins = np.average(res_bins, 1)

            temp_arr = self.readout.readoutDictionary[channel]["Temp [mK]"][-1]
            temp_bins = np.reshape(temp_arr, (self.numBins, int(len(voltage_arr)/self.numBins)))
            temp_bins = np.average(temp_bins, 1)

            self.bin_dictionary[channel]["V [V]"].append(voltage_bins)
            self.bin_dictionary[channel]["R [kohms]"].append(res_bins)
            self.bin_dictionary[channel]["Temp [mK]"].append(temp_bins) 

    def PrintReadout(self, channels_to_print:ArrayLike):
        """
        Prints the voltage(V), resistance(kOhms), and temeprature(K) for the specified channels.

        Parameters
        ----------
        channels_to_print : ArrayLike
            An array specifying the names of the channels to be printed to the console.
            Eg. ["AIN50", "AIN51"]

        Returns
        -------
        None
        """

        print("#################################################################\n")
        print("Stream: {0}".format(self.readout.stream_num))
        for channel in channels_to_print:

            voltage = self.avg_dictionary[channel]["V [V]"][-1]
            temp = self.avg_dictionary[channel]["Temp [mK]"][-1]
            res = self.avg_dictionary[channel]["R [kohms]"][-1]

            print("Channel: {0}; Voltage: {1}V; Resistance: {2}kOhms; Temperature: {3}K".format(channel, str(voltage), str(res), str(temp/1000)))

        print("\n")

    def TestingRoutine(self):
        """
        Specifies the testing routine which is called a the end of every stream in the main streamout class.
        This may be modified along with helper functions to obtain the wanted statistics/data.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.AppendAverageReadout()
        self.AppendBinReadouts()
        self.PrintReadout(self.channels_to_print)

        self.Plot(self.channels_to_plot)

        #Save Data Files
        if self.readout.stream_num % self.saveInterval == 0:
            self.SaveDataFiles() 

    def LoadData(self, channels:ArrayLike, path:str = "./"):
        """
        Loads previous testing data into the `LoadedDataDictionaryAveraged` and `LoadedDataDictionaryBinned` objects.
        The naming convention for testing data must be upheld for this function to work.

        Parameters
        ----------
        channels : ArrayLike
            An array of channels specified by name to be loaded.
            Eg. ["AIN50", "AIN51"]
        
        path : str, default = './'
            The path of the folder which contains the data to be loaded.

        Returns
        -------
        None
        """

        loadDictionary = {"V":"V [V]","R":"R [kohms]","T":"Temp [mK]"}
        
        for channel in channels:
            
            self.LoadedDataDictionaryAveraged[channel] = {"V [V]":[],"R [kohms]":[],"Temp [mK]":[],"Time":[]}
            self.LoadedDataDictionaryBinned[channel] = {"V [V]":[],"R [kohms]":[],"Temp [mK]":[],"Time":[]}
        
            for dataType in loadDictionary: 

                count = 1
                rootFilename = "{0}{1}_{2}".format(path, channel, dataType)
                filename = "{0}_{1}.npz".format(rootFilename, count)
                
                while(os.path.isfile(filename)):
                    
                    raw_data = np.load(filename)

                    self.LoadedDataDictionaryAveraged[channel][loadDictionary[dataType]] = np.append(self.LoadedDataDictionaryAveraged[channel][loadDictionary[dataType]], raw_data["average_" + dataType])
                    self.LoadedDataDictionaryBinned[channel][loadDictionary[dataType]] = np.append(self.LoadedDataDictionaryBinned[channel][loadDictionary[dataType]], raw_data["binned_" + dataType])
                    
                    if dataType == "T": #Appends Time data on last cycle
                        self.LoadedDataDictionaryAveraged[channel]["Time"] = np.append(self.LoadedDataDictionaryAveraged[channel]["Time"], raw_data["average_"+dataType+"_time"])
                        self.LoadedDataDictionaryBinned[channel]["Time"] = np.append(self.LoadedDataDictionaryBinned[channel]["Time"], raw_data["binned_"+dataType+"_time"])

                    count += 1
                    filename = "{0}_{1}.npz".format(rootFilename, count)

        print('Files Loaded.')
                    
    def SaveDataFiles(self):
        """
        Saves the data in `LoadedDataDictionaryAveraged` and `LoadedDataDictionaryBinned` into files that can be later loaded.
        
        Naming Convention: <ChannelName>_<DataType>_<SaveIntervalCounter>
        <ChannelName> - The name of the channel being saved. The channel list is the same as that of the main readout object.
        <DataType> - The data being stored. V: Voltage, R: Resistance, T: Temperature
        <SaveIntervalCounter> - Specifies the the save number of the current instance. Set based on `saveCounter`. Greater the value of `saveInterval`, smaller the number of files for a given time frame.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """     

        self.saveCounter += 1
        split_time_array = self.CalculateSplitTime()
        
        if not os.path.isdir(os.path.realpath(self.savePath)):
            os.makedirs(self.savePath)

        for channel in self.readout.readoutDictionary:
            filename = self.savePath + channel + "_V_" + str(self.saveCounter)
            np.savez(filename, average_V_time = self.avg_dictionary[channel]["Time"], average_V = self.avg_dictionary[channel]["V [V]"], binned_V_time = split_time_array, binned_V = self.bin_dictionary[channel]["V [V]"])
            
            filename = self.savePath + channel + "_R_" + str(self.saveCounter)
            np.savez(filename, average_R_time = self.avg_dictionary[channel]["Time"], average_R = self.avg_dictionary[channel]["R [kohms]"], binned_R_time = split_time_array, binned_R = self.bin_dictionary[channel]["R [kohms]"])
            
            filename = self.savePath + channel + "_T_" + str(self.saveCounter)
            np.savez(filename, average_T_time = self.avg_dictionary[channel]["Time"], average_T = self.avg_dictionary[channel]["Temp [mK]"], binned_T_time = split_time_array, binned_T = self.bin_dictionary[channel]["Temp [mK]"])
         
        self.GenerateDictionaries()

    def RescaleTime(self, time_array:ArrayLike, time_mode:str = "minutes"):
        """
        Rescales time for graphing functions relative to the time of the first data point.

        Parameters
        ----------
        time_array : ArrayLike
            An array containing the time data in Unix time.

        time_mode : str
            String specifying the time interval to be outputted.
            Accepts: "seconds", "minutes", "hours"

        Returns
        -------
        time_array : ArrayLike
            Array with the rescaled time values.
        """  

        if time_mode == "seconds":
            return (time_array - time_array[0])/1
        elif time_mode == "minutes":
            return (time_array - time_array[0])/60
        elif time_mode == "hours":
            return (time_array - time_array[0])/3600
        else:
            return 0

    def GraphData(self, channels:ArrayLike, graph_mode:str = "a", attr:str = "Temp [mK]", time_mode:str = "minutes"):
        """
        Graphs the specified data for the given specified channels. NB: `LoadData()` must be called prior to calling this function.
         
        For more flexibility, the data loaded can be pulled directly from
        `LoadedDataDictionaryAveraged` and `LoadedDataDictionaryBinned`.

        Parameters
        ----------
        channels : Arraylike
            An array of the channels to be read. Each element of the array must be the channel name.
            Eg: ["AIN50", "AIN51", "AIN52"]

        time_mode : str, default = "minutes"
            String specifying the time interval to be outputted.
            Accepts: "seconds", "minutes", "hours"

        attr : str, default = "Temp [mK]"
            String specifying the data type of the data to be plotted.
            Accepts: "V [V]"; "R [kohms]"; "Temp [mK]"

        graph_mode : str, default = "a"
            String specifying the graph to be outputted.
            Modes: 'a' : Average; 's' : Split/Binned; 'as' : "Split and Average"

        Returns
        -------
        None
        """ 

        graphTitle = "null"

        if graph_mode == 'a':
            graphTitle = "Graph Showing Average Data for Specified Channels"
            for channel in channels:
                self.ax.plot(self.RescaleTime(self.LoadedDataDictionaryAveraged[channel]["Time"], time_mode),self.LoadedDataDictionaryAveraged[channel][attr], label = channel)
        
        elif graph_mode == 's':
            graphTitle = "Graph Showing Binned Data for Specified Channels"
            for channel in channels:
                self.ax.scatter(self.RescaleTime(self.LoadedDataDictionaryBinned[channel]["Time"], time_mode),self.LoadedDataDictionaryBinned[channel][attr], label = channel)
        
        elif graph_mode == 'as':
            graphTitle = "Graph Showing Average and Binned Data for Specified Channels"
            for channel in channels:
                self.ax.plot(self.RescaleTime(self.LoadedDataDictionaryAveraged[channel]["Time"], time_mode),self.LoadedDataDictionaryAveraged[channel][attr], label = channel + "_avg")
                self.ax.scatter(self.RescaleTime(self.LoadedDataDictionaryBinned[channel]["Time"], time_mode),self.LoadedDataDictionaryBinned[channel][attr], label = channel + "_bin")
        
        else:
            print("Undefined graph mode.")
            return

        plt.title(graphTitle)
        plt.xlabel("t [" + time_mode + "]")
        plt.ylabel(attr)
        plt.legend()
        plt.show()

#Stream Test
readoutObj = ColdThermometryReadout('56')
readoutObj.Stream()

#Load Data
# testing = Testing()
# testing.LoadData(["AIN56"], './Board Test/')
# testing.GraphData(["AIN56"], "as")