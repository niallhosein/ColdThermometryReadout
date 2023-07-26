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
```
  0.050   6.37650933339e+04  -2.88865353564e+06  -2.26507434131e+00
  0.055   5.21067061211e+04  -1.94906258904e+06  -2.05728687107e+00
  0.060   4.35157796269e+04  -1.49185928308e+06  -2.05699076869e+00
  0.065   3.71238961461e+04  -1.07221325005e+06  -1.87733154352e+00
  0.070   3.26012827893e+04  -7.64638274985e+05  -1.64179672300e+00
  0.075   2.92533099785e+04  -5.92727607619e+05  -1.51964241325e+00
  0.080   2.65636182357e+04  -4.88024137535e+05  -1.46975199901e+00
  0.085   2.43306768135e+04  -4.08840600817e+05  -1.42829775496e+00
  0.090   2.24436802733e+04  -3.48649240536e+05  -1.39809653614e+00
  0.095   2.08203413591e+04  -3.02798241258e+05  -1.38162157975e+00
  0.100   1.94004680946e+04  -2.66199310252e+05  -1.37212828553e+00
  0.110   1.70391291080e+04  -2.09256937987e+05  -1.35090608403e+00
  0.120   1.51654632433e+04  -1.67395884989e+05  -1.32455605717e+00
  0.130   1.36582172879e+04  -1.35468876436e+05  -1.28940355579e+00
  0.140   1.24316125712e+04  -1.10908415477e+05  -1.24900756663e+00
  0.150   1.14216091563e+04  -9.18965804613e+04  -1.20687784712e+00
  0.160   1.05795914102e+04  -7.71234758915e+04  -1.16637360217e+00
  0.170   9.86859570001e+03  -6.55234046094e+04  -1.12872987426e+00
  0.180   9.26110701331e+03  -5.63332832749e+04  -1.09490053132e+00
  0.190   8.73585845713e+03  -4.89710713256e+04  -1.06509321294e+00
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
