# Cold Thermometry Readout

## Structure

## Setup

### Resistor Calibration Curve

### Mapping Calibration Curves for Each Channel
The resistor calibration mapping is specifed in a `.csv` file 

| Left-Aligned  | Center Aligned  | Right Aligned |
| :------------ |:---------------:| -----:|
| col 3 is      | some wordy text | $1600 |
| col 2 is      | centered        |   $12 |
| zebra stripes | are neat        |    $1 |

### Starting a Stream Session
~~~python
readoutObj = ColdThermometryReadout('56')
readoutObj.Stream()
~~~
### Calibration
