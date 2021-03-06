"mirror_power_h.mph" - is the COMSOL file that is used for parametric sweeps of 
the heater power and the convective heat transfer parameter - Figure 6 in the paper.
This file imports data from "input_data.txt"

"mirror_h_function_temperature.mph" - is the COMSOL file used to generate the temperature
transients when the convective heat transfer parameter h is a function of temperature. This file imports
data from "h_data.txt" and "input_data.txt" - It is used to generate Fig. 4 in the paper

"mirror_emissivity.mph" - is the COMSOL file used for parametric sweep of emissivities. This file imports
data from "h_data.txt" and "input_data.txt"

"state_space_export.mph" - is the COMSOL file used to define the model that is used to generate the MATLAB
file "state_space_export.m"

"state_space_export.m" - is the MATLAB file used to generate state-space matrices. These matrices are used in the file 
"model_reduction.m" to perform model order reduction 

"model_reduction.m" - is the MATLAB file used to perform model order reduction. This file load the matrix data 
generated by "state_space_export.m" . This file is used to generate Fig. 9 in the paper

"h_data.txt" - this file is generated by the file "compute_h_varying_temperature.m" (from the folder "preliminary calculations"). This file contains data points
that are used to define dependence of the convective heat transfer coefficient with respect the temperature. This file is imported in COMSOL simulations

"input_data.txt"  - this file is generated by "input_generation.m" (from the folder "preliminary calculations") and imported in COMSOL simulation. This file 
contains the data points defining the heater input (heater on and heater off for predifined time intervals).
