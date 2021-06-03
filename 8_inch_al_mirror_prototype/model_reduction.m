% assumption: 
% both "sssMOR" and "sss" toolboxes are installed and included in MATLAB path 
% the toolboxes can be found here
% https://github.com/MORLab
clear
% this file contrains the matrices exported from COMSOL
load matrices_aluminum_mirror_thick

Am2=sparse(Am);
Bm2=sparse(Bm);
Cm2=sparse(Cm);
Em2=sparse(Em);

[r,n]=size(Cm2)
[~,m]=size(Bm2)

% create the space state-space model
sys = sss(Am2,Bm2,Cm2,[],Em2);

% discretization step
Ts=30; % for best results, the discretization step should be selected on the basis of the step response of the system

%simulation time in discrete-time steps
simulationTime=100

% input for simulation
input2=10*randn(simulationTime,m);

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
model_order=16
[sysr,VW,D]=modalMor(sys,model_order);
 

%simulation of the reduced system
x0_reduced=zeros(model_order,1)
[y,x_,tx]=simBackwardEuler(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,input2,x0_reduced,Ts,Ts,1);

% compare the simulation results
plot(0:Ts:(simulationTime-1)*Ts,298+yo(5,:))
hold on 
plot(0:Ts:(simulationTime-1)*Ts,298+y(5,:),'k')

% blue - original system, red - model order 1, black - model order 8,



% compute the relative simulation error
norm(yo-y,'fro')/norm(yo,'fro')

% results
model_order    =[1,   2,  4, 8, 16, 32, 64, 128 ]
error_modalMor =[0.5658,0.4910,  0.4192, 0.3643, 0.3488, 0.3276, 0.3203, 0.3145]

plot(model_order,error_modalMor,'x-')

disp(sysr)

sys1r=sysr(2,5)
bode(sys1r)

imagesc(abs(inv(sysr.E)*sysr.A))
imagesc(abs(inv(sysr.E)*sysr.B))

eig(inv(sysr.E)*sysr.A)