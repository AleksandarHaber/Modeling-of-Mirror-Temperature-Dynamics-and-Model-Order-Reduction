% This file generates the file 'input_data.txt' that described the heater
% input as a function of time. That is, the heater is "on" for 7179 seconds
% and after it is "off". The file 'input_data.txt' is imported in COMSOL
% simulations to describe the input change
clear, pack, clc

load('total_time_vector')
heater_on_time_instant=find(total_time_vector==7179)
[dim1,~]=size(total_time_vector)
control_vector=zeros(dim1,1)
control_vector(1:heater_on_time_instant,1)=1
control_input_data=[total_time_vector, control_vector]
dlmwrite('input_data.txt',control_input_data)
