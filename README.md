# Cold Thermometry Readout

## Structure

### Settings
~~~python
#Configuable Settings
self.sample_rate:int = 1600
self.scan_amount:int = 10000
self.upper_stage_config:dict = {"UL":100,"LL":1.2,"Bias":5, "CalibratedBias":5}
self.middle_stage_config:dict = {"UL":1.2,"LL":0.15,"Bias":0.2, "CalibratedBias":0.2}
self.lower_stage_config:dict = {"UL":0.15,"LL":0,"Bias":0.02, "CalibratedBias":0.02}
self.BiasSwitchAveragingInterval:int = 30
self.overrideBias:float = 3
self.BiasOutputPort:str = "DAC0"
self.ResistorCalibrationMappingPath:str = "./resistor_calibration_mapping.csv"
self.ReferenceChannel:str = "AIN84"
self.ReferenceChannelResistance:float = 20
self.MemoryBufferSize:int = 3600
    
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
~~~

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
The resistor calibration mapping is specifed in a `.csv` file 

| Channel Name  | Path |
| :------------ | :----|
| AIN48      | ./ResistorCalibrations/AIN48_mean.tbl |
| AIN49      | ./ResistorCalibrations/AIN49_mean.tbl |
| ...        | ... |
| AIN60      | ./ResistorCalibrations/AIN60_mean.tbl |

### Starting a Stream Session
~~~python
readoutObj = ColdThermometryReadout('56')
readoutObj.Stream()
~~~
### Calibration
