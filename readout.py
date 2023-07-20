import os
import sys
import numpy as np
import pandas as pd
from labjack import ljm
from datetime import date, datetime
import plotly.express as px
import time

import matplotlib 
import matplotlib.pyplot as plt

import scipy.integrate as integrate
from scipy import signal
from scipy import stats

import logging

matplotlib.use("tkagg")

logging.basicConfig(filename='coldtherm.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filemode='a')

class ColdThermometryReadout():
    def __init__(self, sample_rate, channels_to_read, scan_amount):
       """
       """
       logging.info("New session started.")

       self.handle = None
       self.setup_names = ''
       self.setup_values = ''
       self.sample_rate = sample_rate
       self.scan_list = None
       self.scan_amount = scan_amount
       self.channel_names = []
       self.readoutDictionary = {}
       self.stream_num = 1

       self.temp_list = []
       self.res_list = []
       self.res_error = []

       #Stage Setup in K and V
       self.upperstage = {"UL":100,"LL":1.2,"Bias":5, "CalibratedBias":5}
       self.middlestage = {"UL":1.2,"LL":0.15,"Bias":0.2, "CalibratedBias":0.2}
       self.lowerstage = {"UL":0.15,"LL":0,"Bias":0.02, "CalibratedBias":0.02}
       self.biasMode = 10 #0: LS, 1: MS, 2: US #3: Override
       self.BiasSwitchAveragingInterval = 30
       
       self.overrideBias = 3
       
       self.BiasOutputPort = "DAC0"

       self.LoadResReadoutGraph()
       self.OpenConnection()
       self.SetupLabJack()
       self.SetChannelsToRead(channels_to_read)
       self.GenerateDictionaries()

       self.TestingModule = Testing(self)

    def BiasSwitch(self, bias:float):
        """
        
        
        """

        ljm.eWriteName(self.handle, self.BiasOutputPort, bias)

    def CheckBiasSwitch(self, ovrr = False):
        """
        
        
        """
      
        if ovrr and self.biasMode == -1:
            return
        elif ovrr:
            self.BiasSwitch(self.overrideBias)
            return

        lowestChannelName = ''
        lowestChannelTemp = 1000
        for channel in self.channel_names:
            
            channeltemp = np.average(self.readoutDictionary[channel]["Temp [mK]"][-self.BiasSwitchAveragingInterval:])/1000 #Temp in K
            if channeltemp < lowestChannelTemp:
                lowestChannelName = channel
                lowestChannelTemp = channeltemp

        if lowestChannelTemp >= self.upperstage["LL"]:
            if self.biasMode != 0:
                self.BiasSwitch(self.upperstage["CalibratedBias"])

        elif lowestChannelTemp < self.upperstage["UL"] and lowestChannelTemp >= self.middlestage["LL"]:
            if self.biasMode != 1:
                self.BiasSwitch(self.middlestage["CalibratedBias"])

        elif lowestChannelTemp < self.lowerstage["UL"]:
            if self.biasMode != 2:
                self.BiasSwitch(self.lowerstage["CalibratedBias"])
       
    def GenerateDictionaries(self):
        """
        ###Testing Function
        This function clears and regenerates all the memory buffers. This is called after the data has been saved to disc.
        """
        self.readoutDictionary.clear()

        for n in self.channel_names:
            self.readoutDictionary[n] = {"V [V]":[],"R [kohms]":[],"Temp [mK]":[],"Time":[]}


    def LoadResReadoutGraph(self, file = os.getcwd() + "\\temp_ref.csv"):
        """
        ###CORE Function
        Loads the Calibration Curve of the thermometer.
        Searches for file in the folder by default. May be specfied.
        """

        file = os.getcwd() + "/temp_ref.csv"
        self.temp_list, self.res_list, self.res_error = np.loadtxt(file, delimiter = ",", skiprows = 0, usecols = (0, 1, 2), unpack = True)


    def OpenConnection(self):
        """
        ###CORE Function
        This function opens a connection to the labjack. Closes all connections to other LabJacks before executing.
        """
        try:
            ljm.closeAll()

            logging.info("Openning connection to LabJack.")

            self.handle = ljm.openS("ANY", "ANY", "ANY")  # Any device, Any connection, Any identifier

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
    
    def SetChannelsToRead(self, channels_to_read):
        """
        ###CORE Function
        Sets the channels to be read.

        Inputs:
        channels_to_read(str) - Input string as '48,56' or 'all' for all channels.
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
    
    def SetupLabJack(self, setup_names=None, setup_values=None):
        """
        # set range, resolution, and settling time for channelsfactor
        # note:
        #   Negative channel: single ended = 199, differential = 1
        #   Range: The instrumentation amplifier in the T7 provides 4 different gains:
        #         x1 (RANGE is ±10 volts), enter 10.0
        #         x10 (RANGE is ±1 volts), enter 1.0
        #         x100 (RANGE is ±0.1 volts), enter 0.1
        #         x1000 (RANGE is ±0.01 volts), enter 0.01
        #   Resolution index = Default (0)
        #   Settling, in microseconds = Auto (0) resource on settling times: https://old3.labjack.com/support/app-notes/SettlingTime
                
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


    def ConvertResistance(self, voltage):
        """"
        Voltage-Resistance Calibraton - Takes an inpt voltage and returns the resistance.

        Input:
        voltage(ArrayLike) - voltage to be converted to a resistance.
        """
        res = voltage
        res = 58.69299438396007*voltage + 0.12049023712211593

        return res

    def Stream(self):
        """
        Periodicity for periodic stream out - 10 correpsonds to print every 10th reading.
        """

        sample_rate = self.sample_rate
        try:
            sample_rate = ljm.eStreamStart(self.handle, int(sample_rate), len(self.channel_names), self.scan_list, sample_rate)
        except:
            ljme = sys.exc_info()[1]
            if ljme.errorCode == 2605:
                print("\n Stream is active. Stopping stream.")
                ljm.eStreamStop(self.handle)
                print("Stream stopped. Restarting Method.")
                self.Stream()

        print("\nStream started with a sample rate of %0.0f Hz." % sample_rate)

            # just a little message
        loop_message = "Press Ctrl+C to stop."
        print("\nStarting %s read loops.%s\n" % (str(self.scan_amount), loop_message))

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
                    temp = np.array(self.ConvertResistance(res))

                    self.readoutDictionary[self.channel_names[k]]["V [V]"].append(v)
                    self.readoutDictionary[self.channel_names[k]]["R [kohms]"].append(res)
                    self.readoutDictionary[self.channel_names[k]]["Temp [mK]"].append(temp)
                    self.readoutDictionary[self.channel_names[k]]["Time"].append(streamtimeutc)
                    pass

                #CHECK FOR BIAS CHANGE HERE
                self.CheckBiasSwitch()

                #CALL TESTING FUNCTIONS HERE
                self.TestingModule.TestingRoutine()

               
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


    def ConvertTemperature(self, resistance, unit="mK"):
        """
        Coonverts resistance to tempermsature.

        unit: 'mK' or 'K'
        """
        factor = 0
        if unit == "mK":
            factor = 1000
        if unit == "K":
            factor = 1
        else:
            pass

        X = factor * np.interp(1000 * resistance, self.res_list, self.temp_list)
        return X

    def CloseConnection(self):
        """
        
        """

        ljm.close(self.handle)
        print("Connection closed.")

    def StreamChannel(self, channel, numStreams):
        """
        Streams the specified channel and returns an array of the streamed averaged values.
        """

        self.SetChannelsToRead(channel)
        voltages = []

        try:
            ljm.eStreamStart(self.handle, int(self.sample_rate), len(self.channel_names), self.scan_list, self.sample_rate)
        except ljm.LJMError:
            ljme = sys.exc_info()[1]
            if ljme.errorCode == 2605:
                #print("\n Stream is active. Stopping stream.")
                ljm.eStreamStop(self.handle)
                #print("Stream stopped. Restarting Method.")
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

        return np.average(voltages) 

class Testing():
    """
    
    """

    def __init__(self, readoutObject:ColdThermometryReadout = None):
        """"
        
        """

        self.readout = readoutObject
        self.avg_dictionary = {}
        self.bin_dictionary = {}

        self.LoadedDataDictionaryAveraged = {}
        self.LoadedDataDictionaryBinned = {}
        
        self.saveInterval = 10
        self.numBins = 10
        self.saveCounter = 0

        self.fig = None
        self.ax = None

        self.GenerateDictionaries()
        self.setup_animation()

    def GenerateDictionaries(self):
        """
        
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
            
            avgVoltage = self.readout.StreamChannel(channel, streamAmount)
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
        ###Testing Function
        """
        #plt.title("Voltage Readout")

        # Updates the animation figure.
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        
    def is_closed(self):

        return not plt.fignum_exists(self.fig.number)
    
    def setup_animation(self):

        #Sets interactive property of the plot.
        plt.ion()

        #Assigns the figure to the game object.
        self.fig = plt.figure()

        #Stores the axes in the plotting object
        self.ax = self.fig.subplots()

    def Plot(self, channels:list, attr:str = "V [V]", graphtitle:str = "Cold Thermometry Test Readout"):
        """
        channels = ["AIN56", "AIN48"]
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
        Calculates the times between splits for the split data array.
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
        Appends the average of a stream to a test dictionary.
        """

        for channel in self.readout.channel_names:
            
            self.avg_dictionary[channel]["V [V]"].append(np.average(self.readout.readoutDictionary[channel]['V [V]'][-1]))
            self.avg_dictionary[channel]["R [kohms]"].append(np.average(self.readout.readoutDictionary[channel]['R [kohms]'][-1]))
            self.avg_dictionary[channel]["Temp [mK]"].append(np.average(self.readout.readoutDictionary[channel]['Temp [mK]'][-1]))
            self.avg_dictionary[channel]["Time"].append(np.average(self.readout.readoutDictionary[channel]['Time'][-1]))
            

    def AppendBinReadouts(self, num_splits:int):
        """
        
        """

        #sr = int(self.sample_rate) ### Converts sample rate to int. Needed to prevent float errors.

        for channel in self.readout.channel_names:
            for j in range(num_splits):
            #a - the lower bound of the split
            #b - the upper bound of the split

                # a = j*int(sr/num_splits)
                # b = ((j + 1) * int(sr/num_splits))
                
                voltage_arr = self.readout.readoutDictionary[channel]["V [V]"][-1]
                voltage_bins = np.reshape(voltage_arr, (num_splits, int(len(voltage_arr)/num_splits)))
                voltage_bins = np.average(voltage_bins, 1)

                res_arr = self.readout.readoutDictionary[channel]["R [kohms]"][-1]
                res_bins = np.reshape(res_arr, (num_splits, int(len(voltage_arr)/num_splits)))
                res_bins = np.average(res_bins, 1)

                temp_arr = self.readout.readoutDictionary[channel]["Temp [mK]"][-1]
                temp_bins = np.reshape(temp_arr, (num_splits, int(len(voltage_arr)/num_splits)))
                temp_bins = np.average(temp_bins, 1)

                self.bin_dictionary[channel]["V [V]"].append(voltage_bins)
                self.bin_dictionary[channel]["R [kohms]"].append(res_bins)
                self.bin_dictionary[channel]["Temp [mK]"].append(temp_bins) 

    def PrintReadout(self, channels_to_print:list):
        """
        
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
        
        """

        ### Caclulate Avergage Data

        self.AppendAverageReadout()
        self.AppendBinReadouts(10)
        self.PrintReadout(['AIN56'])

        self.Plot(["AIN56"])

        #Save Data Files
        if self.readout.stream_num % self.saveInterval == 0:
            self.SaveDataFiles() 

    def LoadData(self, channels:list, path:str = "./"):
        """
        
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

                    t_avg = raw_data["average_"+dataType+"_time"]
                    t_binned = raw_data["binned_"+dataType+"_time"]

                    self.LoadedDataDictionaryAveraged[channel][loadDictionary[dataType]] = np.append(self.LoadedDataDictionaryAveraged[channel][loadDictionary[dataType]], raw_data["average_" + dataType])
                    self.LoadedDataDictionaryBinned[channel][loadDictionary[dataType]] = np.append(self.LoadedDataDictionaryBinned[channel][loadDictionary[dataType]], raw_data["binned_" + dataType])
                    self.LoadedDataDictionaryAveraged[channel]["Time"] = np.append(self.LoadedDataDictionaryAveraged[channel]["Time"], t_avg)
                    self.LoadedDataDictionaryBinned[channel]["Time"] = np.append(self.LoadedDataDictionaryBinned[channel]["Time"], t_binned)

                    count += 1
                    filename = "{0}_{1}.npz".format(rootFilename, count)

        print('Files Loaded.')
                    
    
    def SaveDataFiles(self):
        """"
        ### Core Function
        This Function saves all files to the location of the python file.
        """     

        self.saveCounter += 1
        split_time_array = self.CalculateSplitTime()


        for channel in self.readout.readoutDictionary:
            filename = channel + "_V_" + str(self.saveCounter)
            np.savez(filename, average_V_time = self.avg_dictionary[channel]["Time"], average_V = self.avg_dictionary[channel]["V [V]"], binned_V_time = split_time_array, binned_V = self.bin_dictionary[channel]["V [V]"])
            
            filename = channel + "_R_" + str(self.saveCounter)
            np.savez(filename, average_R_time = self.avg_dictionary[channel]["Time"], average_R = self.avg_dictionary[channel]["R [kohms]"], binned_R_time = split_time_array, binned_R = self.bin_dictionary[channel]["R [kohms]"])
            
            filename = channel + "_T_" + str(self.saveCounter)
            np.savez(filename, average_T_time = self.avg_dictionary[channel]["Time"], average_T = self.avg_dictionary[channel]["Temp [mK]"], binned_T_time = split_time_array, binned_T = self.bin_dictionary[channel]["Temp [mK]"])
         

        self.GenerateDictionaries()

class Analysis():
    """
    
    """

    def __init__(self):
        self.avgdata = np.array([])
        self.avgtime = np.array([])
        self.splitdata = np.array([])
        self.splittime = np.array([])

        self.avgdata2 = np.array([])
        self.avgtime2 = np.array([])
        self.splitdata2 = np.array([])
        self.splittime2 = np.array([])

    def LoadData2(self, root_filename):
        """
        
        """

        count = 1
        filename = "{0}_{1}.npz".format(root_filename, count)

        while(os.path.isfile(filename)):
            raw_data = np.load(filename)

            self.avgdata2 = np.append(self.avgdata2, raw_data["avg"])
            self.avgtime2 = np.append(self.avgtime2, raw_data["avgt"])
            self.splitdata2 = np.append(self.splitdata2, raw_data["split_data"])
            self.splittime2 = np.append(self.splittime2, raw_data["split_time"])

            count += 1
            filename = "{0}_{1}.npz".format(root_filename, count)
        pass

    def LoadData(self, root_filename):
        """
        
        """

        count = 1
        filename = "{0}_{1}.npz".format(root_filename, count)

        while(os.path.isfile(filename)):
            raw_data = np.load(filename)

            self.avgdata = np.append(self.avgdata, raw_data["avg"])
            self.avgtime = np.append(self.avgtime, raw_data["avgt"])
            self.splitdata = np.append(self.splitdata, raw_data["split_data"])
            self.splittime = np.append(self.splittime, raw_data["split_time"])

            count += 1
            filename = "{0}_{1}.npz".format(root_filename, count)
        pass

    def GraphAverageData(self, mode:str):
        """
        mode - time in seconds, minutes or hours
        """

        plt.plot(self.RescaleTime(self.avgtime, mode), self.avgdata)
        plt.plot(self.RescaleTime(self.avgtime2, mode), self.avgdata2)
        plt.plot(self.RescaleTime(self.avgtime, mode), self.avgdata2/self.avgdata, label = "ratio")
        plt.plot(self.RescaleTime(self.avgtime, mode), self.avgdata - self.avgdata2, label = "diff")
        #plt.loglog(abs(np.fft.fft(self.avgdata[1:])), label = "fft data1")

        # import scipy.signal

        # f, pxx = scipy.signal.welch(self.avgdata)
        # plt.plot(f, np.sqrt(pxx))
    
        plt.title("Average Data")
        plt.xlabel("Time/{0}".format(mode))
        plt.ylabel("Voltage/V")
        plt.legend()
        plt.show()

    
    def GraphSplitAverageData(self, mode:str):
        """
         mode - time in seconds, minutes or hours
        """

        plt.plot(self.RescaleTime(self.splittime, mode), self.splitdata, label = "data1")
        plt.plot(self.RescaleTime(self.splittime, mode), self.splitdata2, label = "data2")

        import scipy.signal

        # f, pxx = scipy.signal.welch(self.splitdata, 10)
        # plt.plot(f, pxx)

        plt.title("Split Average Data")
        plt.xlabel("Time/{0}".format(mode))
        plt.ylabel("Voltage/V")
        plt.legend()
        plt.show()

    def GraphAverageDataWithSplitScatter(self, mode:str):
        """
         mode - time in seconds, minutes or hours
        """

        plt.plot(self.RescaleTime(self.avgtime, mode), self.avgdata, color = "r", linewidth = 3, label = "Average Data - data1")
        plt.plot(self.RescaleTime(self.avgtime2, mode), self.avgdata2, color = "r", linewidth = 3, label = "Average Data - data2")
        
        plt.scatter(self.RescaleTime(self.splittime, mode), self.splitdata, label = "Split Data Points - data1")
        plt.scatter(self.RescaleTime(self.splittime2, mode), self.splitdata2, label = "Split Data Points - data2")
        
        plt.title("Average Data with Split Reading Scatter")
        plt.xlabel("Time/{0}".format(mode))
        plt.ylabel("Voltage/V")
        plt.legend()
        plt.show()

    def RescaleTime(self, time_array, mode:str):
        """
         mode - time in seconds, minutes or hours
        """

        if mode == "seconds":
            return (time_array - time_array[0])/1
        elif mode == "minutes":
            return (time_array - time_array[0])/60
        elif mode == "hours":
            return (time_array - time_array[0])/3600
        else:
            return 0
        
    def PowerSpectrum(self, data, dps):
        """
        
        """

        f, pxx = signal.welch(data, dps)
        plt.plot(f, np.sqrt(pxx), label = "Power Spectrum Data 1")
        
        # result = integrate.quad(np.sqrt(pxx), 0, 5)
        # print(result)

        plt.title("Power Spectrum")
        plt.xlabel("Frequency/Hz")
        plt.ylabel("Units")
        plt.legend()
        plt.show()













# #Stream Test
readoutObj = ColdThermometryReadout(1600, '56', 10000)
readoutObj.Stream()






#Load Data
# testing = Testing()
# testing.LoadData(["AIN56", "AIN55"])












#readout.SetScanAmount(3600*24*1.5)
# from threading import Thread, Event
# stop = Event()
# thread = Thread(target=readout.SignalStreamOut, args=(1, 800))
# thread.start()
#readout.Stream()



# analysis = Analysis()
# analysis.LoadData("AIN51__20230707_1600_10")
# analysis.LoadData2("AIN57__20230707_1600_10")
# #analysis.PowerSpectrum(analysis.splitdata, 10)
# analysis.GraphSplitAverageData("hours")

# stop.set(True)



# data = np.load("AIN56__20230705_1600_10_1.npz")

# print(data.files)

# for file in data.files:
#     print("\n\n" + file)
#     print(data[file])
#     print(len(data[file]