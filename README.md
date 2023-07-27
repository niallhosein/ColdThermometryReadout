# Cold Thermometry Readout
This documentation outlines the general steps and important information required to interface with the cold thermometry readout board via software.

## Setup

### Resistor Calibration Curve
The calibration information should be stored in a text file where the respective values in the are separated by spaces as shown below. This is the default format of the `mean.tbl` file. **No header shoulder be present in the file**.

| Temperature(K) | Resistance(Ohms) | dRdT(Ohms/K) | sd |
| :---------: | :---------: |:---------: | :---------:|
|... | ... | ... | ... |

An excerpt of a `mean.tbl` file is shown here:

```
  0.050   6.37650933339D+04  -2.88865353564D+06  -2.26507434131D+00
  0.055   5.21067061211D+04  -1.94906258904D+06  -2.05728687107D+00
  0.060   4.35157796269D+04  -1.49185928308D+06  -2.05699076869D+00
  0.065   3.71238961461D+04  -1.07221325005D+06  -1.87733154352D+00
  0.070   3.26012827893D+04  -7.64638274985D+05  -1.64179672300D+00
```

### Mapping Calibration Curves for Each Channel
It is possible that the resistor connected to each channel of the readout board may have a different temperature-resistance calibration. A mapping file is required to map a calibration curve file to a specific resistor. By default, this file is the `resistor_calibration_mapping.csv` file found in the same directory as the main python script. **A header is expected is in this file.**

Row index 1 is the backup calibration file - the file that will be loaded into a channel in the event that the mapped calibration data cannot be found. The calibration data for other channels can then be explicitly mapped in subsequent rows.

An example calibration mapping is shown below:

| Channel Name  | Path |
| :------------ | :----|
| BACKUP     | ./ResistorCalibrations/BACKUP.tbl |
| AIN48      | ./ResistorCalibrations/AIN48_mean.tbl |
| AIN49      | ./ResistorCalibrations/AIN49_mean.tbl |
| ...        | ... |
| AIN60      | ./ResistorCalibrations/AIN60_mean.tbl |

### Configure the Readout Class
All settings for readout class can be found in the `__init__()` function of the readout class under the heading 'Configurable Settings'. Set the values as desired. A description of each value is given below:

#### Scan Setup
The scan setup is controlled by `sample_rate` and `scan_amount`. `sample_rate` is the number of samples taken by the LabJack per scan per channel. By default it is set to `1600` - a sample frequency of 1600Hz per channel. `scan_amount` is the number of streams to be taken. Set to `-1` for an infinite number of readings.

#### Bias Switching and Stage Configurations
The bias switching functionality is controlled `BiasSwitchAveragingInterval`, `overrideBias`, and `BiasOutputPort`. Sample values are given below:

~~~python
BiasSwitchAveragingInterval:int = 30
self.overrideBias:float = 3 #In volts
self.BiasOutputPort:str = "DAC0"
~~~

There are three stage configurations and their settings are controlled by `upper_stage_config`, `middle_stage_config`, and `lower_stage_config`. Let's consider the middle stage configuration:
    
~~~python
middle_stage_config:dict = {"UL":1.2,"LL":0.15,"Bias":0.2, "CalibratedBias":0.3}

#UL - Upper limit temperature of stage in K
#LL - Lower limit temperature of stage in K
#Bias - The output bias required in volts.
#CalibratedBias - The required DC voltage in volts required to achive the required bias.
~~~

In the above example, the middle stage configuration operates between temperatures of 0.15K and 1.2K, needs to output a bias of 0.2V to best operate the resistors, and requires a DC output of 0.3V from the specified DAC port to achieve the bias of 0.2V.

#### Calibration Mapping
The path to the calibration mapping file must be specified in the `ResistorCalibrationMappingPath` variable. By default, it is set to `"./resistor_calibration_mapping.csv"`.

#### Reference Channel Setup
The reference channel and its fixed resistance value are set in `ReferenceChannel` and `ReferenceChannelResistance`. Sample values are given below:

~~~python
self.ReferenceChannel:str = "AIN84"
self.ReferenceChannelResistance:float = 20 #In kOhms
~~~

#### Miscellaneous
1. `MemoryBufferSize` - The number of scans stored in `readoutDictionary`. After `MemoryBufferSize` scans, `readoutDictionary` is cleared.

#### Calibration
Obtain the calibration equation for the board by calling `CalibrateBoard()` of the `Testing Class`. In the `ConvertResistance()` function of `ColdThermometryReadout`, set the calibration equation obtained as the output. A sample is given below:

~~~python
res = 58.69299438396007 * voltage + 0.12049023712211593
return res
~~~

*This should be a linear equation relating a voltage to a resistance.*

## Testing
The `Testing` class contains all of the functions that are not necessary for the board to operate. This includes: live readout(graphical and to console), statistics functions, testing routines, and data analysis variables. Of importance, it also contains the `CalibrateBoard()` function.

### Setup
All settings for the `Testing` class are found in the `__init__()` function under the section 'Configurable Settings'. The categories are given below:

#### Live Print and Plotting Channels
The channels to be printed and plotted are controlled by `channels_to_print` and `channels_to_plot` respectively, each of which is a list of the channel names for the respective operation. A sample is given below:

~~~python
channels_to_print = ["AIN40", "AIN41", "AIN42"]
channels_to_plot = ["AIN56"]
~~~

#### Save Settings
The save settings are controlled by `saveInterval` and `savePath`. `saveInterval` gives the number of streams which must occur before data is saved to disc. `savePath` gives the path of the folder where the files will be saved.

#### Miscellaneous
1. `numBins` - This is the number of bins to create for each scan. By default, this is set to `10` and with a `scan_rate` of `1600`, this results in 10 bins of size 160. 

### Testing Routine
After each scan, the `TestingRoutine()` function is called. This function calls the necessary functions required for live streaming and saving data for later analysis. A sample testing routine is given below:

~~~python
self.AppendAverageReadout()
self.AppendBinReadouts()
self.PrintReadout(self.channels_to_print)

self.Plot(self.channels_to_plot)

#Save Data Files
if self.readout.stream_num % self.saveInterval == 0:
    self.SaveDataFiles()
~~~

### Calibrate Board
The `CalibrateBoard()` function automatically calibrates the board and outputs a linear equation relating an output voltage to a resistance value. This function prompts the user to set the input resistance of the board to fixed values of their choice, and then streams for a specified number of scans. Using the data obtained from these scans, and the input resistances, it performs linear regression analysis on the data to obtain the calibration equation. 

### Loading Data
Data is saved as `.npz` files which allows multiple `ArrayLike` objects to be stored with ease in a single file. This data can be loaded at a later time by calling the `LoadData` function of the `Testing` class.

~~~python
testing = Testing()
testing.LoadData(["AIN56", "AIN40", "AIN41"], './Data Directory/')
~~~

The loaded data is stored in two dictionaries, `LoadedDataDictionaryAveraged` and `LoadedDataDictionaryBinned`, which store the averaged and binned data respectively.

### Graphing Data
The `GraphData` function of the `Testing` class allows for quick graphing of data. This function has three modes which are set by the `graph_mode` parameter:

1. `a` - Graphs the averaged data for the specified channels as a line plot
2. `s` - Grapahs the binned data for the specified channels as a scatter plot
3. `as` - Graphs the averaged data for the specified channels as a line plot and the binned data for the specified channels as a scatter plot.

**NOTE: `LoadData()` must be called first before using this function.**

*Alternatively, and for more flexibility, the saved data can be loaded, and the raw data required can be obtained from the `LoadedDataDictionaryAveraged` and `LoadedDataDictionaryBinned` dictionaries.*

## Usage
### Starting a Streaming Session
A stream session can be started by creating an instance of the `ColdThermometryReadout` class and calling `Stream()`. **Note: The channels to be read is required when creating an instance of `ColdThermometryReadout`.**

~~~python
readoutObj = ColdThermometryReadout('40,41,42,43') #Use 'all' for all channels.
readoutObj.Stream()
~~~
