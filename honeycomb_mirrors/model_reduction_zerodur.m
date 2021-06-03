% assumption: 
% both "sssMOR" and "sss" toolboxes are installed and included in MATLAB path 
% the toolboxes can be found here
% https://github.com/MORLab
clear
% this file contrains the matrices exported from COMSOL
load matrices_zerodur

Am2=sparse(Am);
Bm2=sparse(Bm);
Cm2=sparse(Cm);
Em2=sparse(Em);

[r,n]=size(Cm2)
[~,m]=size(Bm2)

% create the space state-space model
sys = sss(Am2,Bm2,Cm2,[],Em2);

% discretization step
Ts=4; % for best results, the discretization step should be selected on the basis of the step response of the system

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
model_order=3200
[sysr,VW,D]=modalMor(sys,model_order);
%simulation of the reduced system
x0_reduced=zeros(model_order,1)
[y,x_,tx]=simBackwardEuler(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,input2,x0_reduced,Ts,Ts,1);


% compare the simulation results
plot(0:Ts:(simulationTime-1)*Ts,yo(5,:))
hold on 
plot(0:Ts:(simulationTime-1)*Ts,y(5,:),'g')

% blue - original system, red - model order 10, black - model order 100,
% magenta - model order 200, green - model order 400


% compute the relative simulation error
norm(yo-y,'fro')/norm(yo,'fro')

% results
figure(1)
hold on
model_order    =[10,     25,         50,       100,       200,     400, 800, 1600, 3200]
error_modalMor =[0.7619, 0.7425,   0.7182,  0.1128,  0.0850, 0.0411, 0.0251,  0.0149,0.0038 ]

plot(model_order,error_modalMor,'mx-')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results
figure(1)
hold on
model_order    =[10,     25,         50,       100,       200,     400, 800, 1600, 3200]
error_modalMor =[0.7967,0.7831, 0.7628,0.1206, 0.0945, 0.0438,  0.0274, 0.0160, 0.0043]
plot(model_order,error_modalMor,'kx-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




sys1r=sysr(17,2)
bode(sys1r)



imagesc(abs(inv(sysr.E)*sysr.A))


