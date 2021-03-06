Explanation of the files in this folder:

"temperature_conversion_final.m" - is used to convert raw measurements from the file "Thermistor_Data_21-04-02_1033.lvm"
into temperature reading. This files implements the least squares method and inverts the voltage divider equation. 

"input_generation.m" - is used to generate the file "input_data.txt" that stores the data points of the time function
describing the heater input in the experiment (heater on and heater off for predifined time intervals).

"total_time_vector.mat" - is used by the file "input_generation.m". This is a MATLAB workspace file that contains
the data to define the input time function in "input_data.txt"

"input_data.txt"  - this file is generated by "input_generation.m" and imported in COMSOL simulation. This file 
contains the data points defining the heater input (heater on and heater off for predifined time intervals).

"compute_h_fixed_temperature.m" - this file is used to compute the constant value of the convective heat transfer parameter h

"compute_h_varying_temperature.m" - this file is used to compute the convective heat transfer parameter h as a function of temperature 

"h_data.txt" - this file is generated by the file "compute_h_varying_temperature.m" . This file contains data points
that are used to define dependence of the convective heat transfer coefficient with respect the temperature. This file is imported
in COMSOL simulations




