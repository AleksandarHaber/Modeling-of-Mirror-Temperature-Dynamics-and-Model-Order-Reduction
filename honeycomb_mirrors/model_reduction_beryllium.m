% assumption: 
% both "sssMOR" and "sss" toolboxes are installed and included in MATLAB path 
% the toolboxes can be found here
% https://github.com/MORLab
clear
% this file contrains the matrices exported from COMSOL
load matrices_beryllium.mat

Am2=sparse(Am);
Bm2=sparse(Bm);
Cm2=sparse(Cm);
Em2=sparse(Em);

[r,n]=size(Cm2)
[~,m]=size(Bm2)

% create the space state-space model
sys = sss(Am2,Bm2,Cm2,[],Em2);

% discretization step
Ts=10; % for best results, the discretization step should be selected on the basis of the step response of the system

%simulation time in discrete-time steps
simulationTime=100

% input for simulation
input2=0.5*rand(simulationTime,m);

% initial state in the reduced coordinates
% since the initial condition is zero, we can select zero initial
% conditions for both reduced and large-scale systems, otherwise, we would need 
% to do the following transformation: xo_large_scale=VW*x0_reduced, where VW
% is a transformation matrix coming from the model reduction algorithm

x0_large_scale=zeros(n,1);

% simulation of the full-size system
[yo,xo_,txo]=simBackwardEuler(Am2,Bm2,Cm2,sparse(r,m),Em2,input2,x0_large_scale,Ts,Ts,1);

% model reduction 
% this is the model order for reduction
model_order=10
[sysr,VW,D]=modalMor(sys,model_order);
%simulation of the reduced system
x0_reduced=zeros(model_order,1)
[y,x_,tx]=simBackwardEuler(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,input2,x0_reduced,Ts,Ts,1);


% compare the simulation results
plot(0:Ts:(simulationTime-1)*Ts,yo(5,:))
hold on 
plot(0:Ts:(simulationTime-1)*Ts,y(5,:),'k')

% blue - original system, red - model order 1, black - model order 10,



% compute the relative simulation error
norm(yo-y,'fro')/norm(yo,'fro')

% results
model_order    =[1,     2,         3,       4,       5,     10, 20, 40, 80, 160]
error_modalMor =[0.0457,   0.0451,  0.0447, 0.0446, 0.0444,  0.0256,  0.0199,0.0184, 0.0023,0.0021]

% results - erased 2,3,4 to empty the space
figure(1)
hold on
model_order    =[1,           5,     10, 20, 40, 80, 160]
error_modalMor =[0.0457,  0.0444,  0.0256,  0.0199,0.0184, 0.0023,0.0021]


plot(model_order,error_modalMor,'rx-')



sys1r=sysr(17,2)
bode(sys1r)


