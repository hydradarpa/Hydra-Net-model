%% INTEGRATION PARAMETERs
clear;       % remove previous varibles
close all;

T=300;        % total time, mS
dt=0.01;     % time step, ms

Df=1/dt;            % Delta function approximation

Tframe=0;    % movie
dTframe=50; % ms
frame=10;

figure('units','normalized','outerposition',[0 0 0.8 0.8]); % show figure window
%%

%% NETWORK PARAMETERS
N=441;        % Number of neurons
%%

%% CONNECTIVITY
gEE_mean=7;  % mA/cm^2 (current synapses)
%gEE_var=0.001;   % variance for synaptic connections
tau_EE=5.4;      % ms, AMPA current 5.4

S_EE=M_tor(sqrt(N))*gEE_mean;
%%

%% NEURON PARAMETERS
C=1;        % muF/cm^2
VNa=50;     % mV
gNa=120;    % mS/mm^2
gK=36;      % mS/mm^2
VK=-77;     % mV
gL=0.3;     % mS/mm^2
VL=-54.4;   % mV

Ca_base=0.1; % mM, baseline Ca concentration
Ca_jump=0.5; % mM, jump of Ca after spike
tau1_Ca=100; % Ca rise time
tau2_Ca=200; % Ca decay time

Temp=18;  % Temperature in C, standard preparation
phi=3^((Temp-6.3)/10);

V_spike=0;            % mV, threshold for spike
t_ref=5;              % refractory period for spike recordings
%%

%% ICs
%Constant ICs
%
V(1:N)=-65.15;   % mV
m(1:N)=0.05;     % 1
n(1:N)=0.31;     % 1
h(1:N)=0.60;     % 1

Ca(1:N)=Ca_base;       % baseline Ca concentration
dCa(1:N)=0;
Ca_matrix(round(T/dt),1:N)=0;   % changes in Ca in all neurons
Ca_repr=Ca_base;

%}
% Take the last ICs from single HH simulation
%{
V(1:N)=V(end);   % mV
m(1:N)=m(end);     % 1
n(1:N)=n(end);     % 1
h(1:N)=h(end);     % 1
h(N+1:end)=[];   
V(N+1:end)=[];
m(N+1:end)=[];
n(N+1:end)=[];
%}

%REPRESENTATIVE CELLS
repr=round(N/2);             % number of representative neuron
V_repr=zeros(1);    % mV

% SYNAPTIC VARIABLEs
I_EE=zeros(N,1);

% FIRINGS
firings=[];           % spike timings
fired=[];             % indexes of spikes fired at t
V_sp=10*ones(N,1);   % vector of times of elements that did spike
fired_delta(1:N)=0;   % vector of fired delta-function
%%

%% STIMULATION
IN(1:N)=0;        %mA/cm^2
% stimulated neuron
IN(repr)=8;        %mA/cm^2

ts=5;           % end of stimulation, ms
%%

%% TIME INTEGRATION LOOP
for t=1:1:round(T/dt)
    
% Auxilary functions
%
alfa_m=0.1*(V+40)./(1-exp(-0.1*(V+40)));
beta_m=4.*exp(-0.0556*(V+65));

alfa_n=0.01*(V+55)./(1-exp(-0.1*(V+55)));
beta_n=0.125*exp(-0.0125*(V+65));

% shifted
alfa_h=0.07*exp(-0.05*(V+65));
beta_h=1./(1+exp(-0.1*(V+35)));

% Gating variables
dmdt=(alfa_m.*(1-m)-beta_m.*m)*phi;
m=dmdt*dt + m;
dndt=(alfa_n.*(1-n)-beta_n.*n)*phi;
n=dndt*dt + n;
dhdt=(alfa_h.*(1-h)-beta_h.*h)*phi;
h=dhdt*dt + h;
%}

% Voltage
dVdt=1/C*(-gL.*(V-VL) -gNa.*m.^3.*h.*(V-VNa)-gK*n.^4.*(V-VK) +IN +I_EE');
V=dVdt*dt + V;

% Ca neural dynamics
Ca=dCa*dt + Ca;
dCa=(dt/tau1_Ca/tau2_Ca).*( Ca_jump*fired_delta./K(1/tau1_Ca,1/tau2_Ca)-(tau1_Ca +tau2_Ca).*dCa -(Ca-Ca_base) ) + dCa;
Ca_matrix(t,:)=Ca;         % Ca matrix
Ca_repr(t)=Ca(repr);


% Representative cell
V_repr(t)=V(repr);
I_repr_ext(t)=IN(repr);
I_repr_EE(t)=I_EE(repr);

% Synaptic variable
I_EE=(-I_EE/tau_EE + Df.*sum(S_EE(:,fired),2) )*dt + I_EE;

% Stimulation

if t*dt>ts
    IN(repr)=0;
end

% Firings processing

% Debug
if isempty(find(V>V_spike))==0
    a=1;
end

fired=[];
fired=find(V>V_spike);

INT_E=(intersect(find(V_sp<t_ref),fired))';   % intersection V_sp<2 and fired_E
fired=setxor(fired,INT_E);                    % remove fired_E elements that intersect with VSOMA_sp<2

V_sp=V_sp + dt;                               % update the t* vector
V_sp(fired)=0;
            
t_fired(1:length(fired))=t;
        
% check the format!!!
if isrow(fired)==0
     fired=fired';
end
    
if isrow(t_fired)==0
     t_fired=t_fired';
end
    
spikes=horzcat(t_fired',fired');   
firings=vertcat(firings,spikes);

t_fired=[];
spikes=[];

fired_delta(:)=0;           % vector of current delta functions of all neurons
fired_delta(fired)=Df;      

%}
%end
%%

%% FINAL PLOT / ONLINE PLOT

if (t-1)*dt==Tframe
    
subplot(2,2,1);
if isempty(firings)==0
plot(firings(:,1)*dt,firings(:,2),'.','MarkerSize',15);
end
ylabel('Cell index');
ylim([0 N]);
set(gca,'FontSize',20);             % set the axis with big font
title(sprintf('HH population, T=%d ms',round(Tframe)));
set(gca,'FontSize',20);             % set the axis with big font
box off;

subplot(2,2,2);
plot((1:t)*dt,V_repr);
ylabel('Voltage (mV)');
set(gca,'FontSize',20);             % set the axis with big font
title('Representative neuron');
set(gca,'FontSize',20);             % set the axis with big font
box off;

subplot(2,2,3);
imagesc((reshape(V,sqrt(N),sqrt(N)))',[-80 40]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
title('Voltage distribution (mV)');
colorbar;
box off;

subplot(2,2,4);
imagesc((reshape((IN +I_EE'),sqrt(N),sqrt(N)))',[0 30]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
colorbar;
title('Synaptic input (\muA/cm^2)');
box off;

MOV(frame)=getframe(gcf);

frame=frame+1;          % counter for the movie
Tframe=Tframe + dTframe;
%%

end

end

% Save the movie
%movie2avi(MOV,'HH_net.avi','fps',10,'quality',1);

%%