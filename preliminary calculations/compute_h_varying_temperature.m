% This code is used to compute the free-convection heat transfer
% coefficient. Here we assume that the surface temperature is changing
% Author: Aleksandar Haber, May 06, 2021
% This MATLAB file generates the file 'h_data.txt' that stored data points
% representing the dependence of the convective heat transfer parameter on
% the surface temperature. The file 'h_data.txt' is imported in COMSOL
% simulations
clear, pack, clc

% initial temperature
T0=298 % THIS VALUE YOU CAN ADJUST USING THE EXPERIMENTAL DATA

% gravitational acceleration constant
g=9.81

% length of the vertical plate 
L=0.2032

% data for interpolation in celsius
% kinematic viscosity - nu
% thermal conductivity - k
% Prandtl number - Pr
% Table A.6-SI, page 556 from "Heat Transfer", Fourth Edition, by Alan J. Chapman
temperature=[0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200]

nu=[13.31 14.19 15.09 16.01 16.96 17.92 18.90 19.90 20.92 21.96 23.02 24.10 25.19 26.31 27.44 28.58 29.75 30.93 32.13 33.34 34.57]*10^(-6)
k=[24.08 24.87 25.64 26.38 27.10 27.81 28.52 29.22 29.91 30.59 31.27 31.94 32.61 33.28 33.94 34.59 35.25 35.89 36.54 37.18 37.81]*10^(-3)
Pr=[0.718 0.716 0.713 0.712 0.710 0.709 0.708 0.707 0.706 0.705 0.704 0.704 0.703 0.702 0.702 0.701 0.701 0.700 0.700 0.699 0.699]


for i=1:200
% maximal temperature
% temperature of the surface
Tmax=T0+i

% mean film temperature
Tf=(Tmax+T0)/2

% temperature difference
deltaT=Tmax-T0

% mean film temperature in celsius
Tfcelsius=Tf-273.15

nu_int = interp1(temperature,nu,Tfcelsius)
k_int  = interp1(temperature,k,Tfcelsius)
Pr_int = interp1(temperature,Pr,Tfcelsius)

% volumetric thermal expansion coefficient 
beta = 1/T0

%Grashof number
Gr_int=(g*beta*deltaT*L^3)/(nu_int^2)

%Rayleigh number
Ra_int=Pr_int*Gr_int

%Nusselt number
Nu_int(i)=0.68+(0.670*(Ra_int)^(1/4))/((1+(0.492/Pr_int)^(9/16))^(4/9))

%average free-convection heat transfer coefficient
h_loop(i)=Nu_int(i)*k_int/L
 
end
temperature_loop=1:200;
temperature_loop=T0+temperature_loop;

figure(1)
plot(temperature_loop,h_loop)

temperature_loop=[288 temperature_loop]
h_loop=[h_loop(1) h_loop]


h_data=[temperature_loop', h_loop']

dlmwrite('h_data.txt',h_data)


