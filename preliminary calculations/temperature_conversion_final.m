% This file is used to convert the raw measurements given in the file:
% "Thermistor_Data_21-04-02_1033.lvm" into temperature readings

%It is Temperature Vs Resistance curve.
clear,pack,clc

% data taken from 
% https://media.digikey.com/pdf/Data%20Sheets/Ametherm%20PDFs/ACC-101_Dwg.pdf

R = [75780, 39860, 21860, 12460, 7352.8, 4481.5, 2812.8, 2252, 1814.4, 1199.6, 811.40, 560.30, 394.55, 282.63, 206.13, 152.75, 114.92, 87.671, 67.77, 52.983, 41.881];
T = 273.15+[-40,-30,-20,-10,0,10,20,25,30,40,50,60,70,80,90,100,110,120,130,140,150];

T=T(6:end)
R=R(6:end)

%plot(R,T)

% getting the calibration constants:a,b,c
% 1/T=a+b*log(R)+c*(log(R))^3
% y=a+b*x+c*x^3
% x=log(R)
% y=1/T
% Steinhart–Hart equation
%https://en.wikipedia.org/wiki/Thermistor
% we use the least-squares method

x=log(R)
y=1./T
s1=numel(R)

for i=1:s1
    A(i,1)=1
    A(i,2)=x(i)
    A(i,3)=x(i)^3
end

parameters=inv(A'*A)*A'*y'

data1 = csvread('Thermistor_Data_21-04-02_1033.lvm');

[s2,~]=size(data1(:,1))
measurement1=data1(:,4);
time_vector=data1(:,1);

% get the temperature readings

for i=1:s2
   
   measurement2(i)=2000*((4.98/(4.98-measurement1(i)))-1);
   T(i)=1/(parameters(1)+parameters(2)*log(measurement2(i))+parameters(3)*(log(measurement2(i)))^3); 

end

figure(1)
hold on
plot(time_vector,T,'b')




















time_vector2=0:30:2790
%temp_comsol=comsol_data(4,3:end)
plot(time_vector2,comsol_data,'-')







% V^2/R = 12V^2/37ohm = 144/37 = 3.9 W
% power computation

resistance=37;
((12)^2)/resistance
((0.8*12)^2)/resistance
((0.6*12)^2)/resistance
((0.4*12)^2)/resistance
((0.2*12)^2)/resistance


voltage=[0.2*12 0.4*12 0.6*12 0.8*12 12]
temperature=[296.7403 297.9817 299.7257 302.3536 306.4667]

plot(voltage, temperature)



figure(5)
plot(data1(:,1),data1(:,9))



