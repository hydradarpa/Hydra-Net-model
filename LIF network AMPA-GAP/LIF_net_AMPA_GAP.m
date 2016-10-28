%% INTEGRATION PARAMETERS
clear;       % remove previous varibles
%close all;

T=110;        % total time, mS
dt=0.01;     % time step, ms

%time=(1:1:(T/dt)).*dt;      % TIME VECTOR
Df=1/dt;            % Delta function approximation

Tframe=0;    % movie
dTframe=10; % ms
frame=10;

figure('units','normalized','outerposition',[0 0 0.8 0.8]); % show figure window
%%

%% NETWORK PARAMETERS
N=676;        % Number of neurons
%%

%% CONNECTIVITY
% Synapses
gEE_mean=20;  % mA/cm^2 (current synapses)
tau_EE=5.4;      % ms, AMPA current 5.4
S_EE=M_tor(sqrt(N))*gEE_mean;

% GAP junctions, the same location as synapses
g_GAP=0.03;
l=zeros(N,4);
for i=1:1:N 
   l(i,1:length(find(S_EE(i,1:end)>0)))=find(S_EE(i,1:end)>0);
end


% Representative cell
repr=round(N/2);             % number of representative neuron

%%

%% NEURON PARAMETERS
taum=10;        % muF/cm^2
R=1;            % K Ohm
VL=-60;         % mV
Vreset=-60;     % -6
VT=-55;         % mV, threshold for spike
t_ref=15;        % refractory period

% Ca dynamics
Ca_base=0.1; % mM, baseline Ca concentration
Ca_jump=0.5; % mM, jump of Ca after spike
tau1_Ca=100; % Ca rise time
tau2_Ca=200; % Ca decay time


%%

%% STIMULATION
IN(1:N)=0;          % mA/cm^2
% stimulated neuron
IN(repr)=15;        % mA/cm^2
ts=5;              % end of stimulation, ms
%%

%% ICs
%Constant ICs
V(1:N)=VL;   % mV
V_repr=zeros(1); % mV

% Ca dynamics
Ca(1:N)=Ca_base;            % current Ca concentration
dCa(1:N)=0;
Ca_matrix(round(T/dt),1:N)=0;

% SYNAPTIC VARIABLES
I_EE=zeros(N,1);                % synaptic current
I_GAP=zeros(N,1);               % GAP-junction current

% FIRINGS
firings=[];           % spike timings
fired=[];             % indexes of spikes fired at t
V_sp=50*ones(N,1);    % vector of times of elements that did spike

fired_delta(1:N)=0;   % vector of fired delta-function
%%

%% TIME INTEGRATION LOOP
for t=1:1:round(T/dt)

% Stimulation
if t*dt>ts
    IN(repr)=0;
end
    
% Connections
I_EE=(-I_EE/tau_EE + Df.*sum(S_EE(:,fired),2) )*dt + I_EE;      % Synaptic current
I_GAP=g_GAP*(sum(V(l),2)-4*V');                                 % GAP-junction current (generalize for open borders!)

% Voltage
dVdt=1/taum*(-R.*(V-VL) +IN +I_EE' +I_GAP');
V=dVdt*dt + V;


% Ca neural dynamics
Ca=dCa*dt + Ca;
dCa=(dt/tau1_Ca/tau2_Ca).*( Ca_jump*fired_delta./K(1/tau1_Ca,1/tau2_Ca)-(tau1_Ca +tau2_Ca).*dCa -(Ca-Ca_base) ) + dCa;
Ca_matrix(t,:)=Ca;         % Ca matrix


% Representative cell
V_repr(t)=V(repr);
I_repr_ext(t)=IN(repr);
I_repr_EE(t)=I_EE(repr)+I_GAP(repr);

if V_repr(t)>=VT
    V_repr(t-1)=10;     % add stick to the repr neuron
end

% FIRINGS
fired=[];
fired=find(V>=VT);
V(fired)=Vreset;

% refractory period
V_sp=V_sp + dt;                               
V_sp(fired)=0;
V(find(V_sp<t_ref))=Vreset;          

t_fired(1:length(fired))=t;

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
title(sprintf('LIF population, T=%d ms',round(Tframe)));
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
imagesc((reshape(V,sqrt(N),sqrt(N)))',[Vreset VT]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
title('Voltage distribution (mV)');
colorbar;
box off;

subplot(2,2,4);
imagesc((reshape((IN +I_EE' +I_GAP'),sqrt(N),sqrt(N)))',[0 gEE_mean*4]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
colorbar;
title('Synaptic and GAP-junction input (\muA/cm^2)');
box off;

MOV(frame)=getframe(gcf);

frame=frame+1;          % counter for the movie
Tframe=Tframe + dTframe;
%%

end

end

% Save the movie
%movie2avi(MOV,'LIF_net.avi','fps',10,'quality',1);

%% Plot delay distribution
%{
figure('units','normalized','outerposition',[0 0 0.7 0.7]); % show figure window

histogram((firings(:,1)*dt-firings(1,1)*dt),10) % histogram with respect to the 1st spike

set(gca,'Fontsize',30)
xlabel('Spike delay (ms)')
ylabel('Count')
title('LIF network delay distribution (first spike)')
%}

%%