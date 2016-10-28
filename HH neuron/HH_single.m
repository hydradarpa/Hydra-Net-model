%% NUMERICAL PARAMETERS
clear;       % reset everything before run

T=5;        % total time, mS
dt=0.01;     % time step, ms
%%

%% MODEL PARAMETERS
C=1;        % muF/cm^2
VNa=50;     % mV
gNa=120;    % mS/mm^2
gK=36;      % mS/mm^2
VK=-77;     % mV
gL=0.3;     % mS/mm^2
VL=-54.4;     % mV

Temp=18;  % Temperature in C, standard preparation
phi=3^((Temp-6.3)/10);

I=8;        %mA/cm^2
%%

%%
%% INITIAL CONDITIONS (rest state)
V(1:1:round(T/dt))=-65.15;   % mV
m(1:1:round(T/dt))=0.05;     % 1
n(1:1:round(T/dt))=0.31;     % 1
h(1:1:round(T/dt))=0.60;     % 1
%%


%% SIMULATION
t=(0:1:round(T/dt))*dt; % time vector

for i=1:1:round(T/dt)
    
% Auxilary functions
alfa_m=0.1*(V(i)+40)/(1-exp(-0.1*(V(i)+40)));
beta_m=4*exp(-0.0556*(V(i)+65));

alfa_n=0.01*(V(i)+55)/(1-exp(-0.1*(V(i)+55)));
beta_n=0.125*exp(-0.0125*(V(i)+65));
alfa_h=0.07*exp(-0.05*(V(i)+65));
beta_h=1./(1+exp(-0.1*(V(i)+35)));

% Gating variables
dmdt=(alfa_m*(1-m(i))-beta_m*m(i))*phi;
m(i+1)=dmdt*dt + m(i);
dndt=(alfa_n*(1-n(i))-beta_n*n(i))*phi;
n(i+1)=dndt*dt + n(i);
dhdt=(alfa_h*(1-h(i))-beta_h*h(i))*phi;
h(i+1)=dhdt*dt + h(i);

% Voltage
dVdt=1/C*(-gNa*m(i)^3*h(i)*(V(i)-VNa)-gK*n(i)^4*(V(i)-VK)-gL*(V(i)-VL)+I);
V(i+1)=dVdt*dt + V(i);
    
end
%%

%% PLOT
figure;
plot(t,V);
set(gca,'FontSize',30);             % set the axis with big font
xlabel('time (ms)');
ylabel('Voltage (mV)');
%%

