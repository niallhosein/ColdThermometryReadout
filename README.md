# Cold Thermometry Readout

## Structure

### Settings

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

#### sample_rate : int
This is the number of samples taken by the LabJack per scan per channel. By default it is set to `1600` - a sample frequency of 1600Hz per channel.

#### scan_amount : int
The number of streams to be taken. Set to `-1` for an infintie number of readings.

#### Bias Switching and Stage Configurations
The bias switching functionality is controlled `BiasSwitchAveragingInterval`, `overrideBias`, and `BiasOutputPort`. Sample values are given below:

~~~python
BiasSwitchAveragingInterval:int = 30
self.overrideBias:float = 3 #In volts
self.BiasOutputPort:str = "DAC0"
~~~

There are three stage configurations and there settings are controlled by `upper_stage_config`, `middle_stage_config`, and `lower_stage_config`. Let's consider the middle stage configuration:
    
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

### Starting a Stream Session
~~~python
readoutObj = ColdThermometryReadout('56')
readoutObj.Stream()
~~~
### Calibration
