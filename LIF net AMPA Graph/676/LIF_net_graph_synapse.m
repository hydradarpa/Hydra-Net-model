%% INTEGRATION PARAMETERS
clear;       % remove previous varibles
close all;

T=1000;        % total time, mS
dt=0.01;     % time step, ms

%time=(1:1:(T/dt)).*dt;      % TIME VECTOR
Df=1/dt;            % Delta function approximation

Tframe=0;    % movie
dTframe=100; % ms
frame=100;

figure('units','normalized','outerposition',[0 0 0.8 0.8]); % show figure window
%%

%% NETWORK PARAMETERS
N=676;        % Number of neurons
%%

%% CONNECTIVITY
gEE_mean=10;  % mA/cm^2 (current synapses)  # 37
tau_EE=5;      % ms, AMPA 5.4

load('A_676.mat')

S_EE=A*gEE_mean;    % connectivity matrix

% Representative cell
repr=round(N/2);             % number of representative neuron

%%

%% NEURON PARAMETERS
taum=10;        % muF/cm^2
R=1;            % K Ohm
VL=-60;         % mV
Vreset=-60;     % -6
VT=-55;         % mV, threshold for spike
t_ref=10;        % refractory period, ms

%%

%% STIMULATION
IN(1:N)=0;          % mA/cm^2

% stimulated neuron
%stim_1=randperm(N,10);
%stim_1=(350:375);          % (325:375)
%stim_1=datasample(find(sum(A)==1),7);    % pick up receptor neurons
stim_1=339;                % central neuron

stim_2=(320:345);          % critical
stim_3=randperm(N,10);
stim_4=randperm(N,10);

IN(stim_1)=0;        % mA/cm^2
IN(stim_2)=0;        % mA/cm^2
IN(stim_3)=0;        % mA/cm^2
IN(stim_4)=0;        % mA/cm^2

ts1_start=1;         % start of first stimuli
ts1_end=51;          % end of first stimuli

ts2_start=1000;       % start of second stimuli
ts2_end=1100;         % end of second stimuli

ts3_start=500;       % start of third stimuli
ts3_end=520;         % end of third stimuli

ts4_start=750;       % start of third stimuli
ts4_end=770;         % end of third stimuli

%}
%%

%% ICs
%Constant ICs
V(1:N)=VL;   % mV
%V(1:N)=random('Normal',VL,10,1,N);

V_repr=zeros(1); % mV

% SYNAPTIC VARIABLES
I_EE=zeros(N,1);

% FIRINGS
firings=[];           % spike timings
fired=[];             % indexes of spikes fired at t
V_sp=50*ones(N,1);    % vector of times of elements that did not spike

fired_delta(1:N)=0;   % vector of fired delta-function
%%

%% TIME INTEGRATION LOOP
for t=1:1:round(T/dt)

% Stimulation

if t*dt>ts1_start
    if t*dt<ts1_end
       IN(stim_1)=10;
    else
        IN(stim_1)=0;
    end
end

IN_stim(t)=IN(stim_1(1));

%{

if t*dt>ts2_start
   if t*dt<ts2_end
       IN(stim_2)=10;
   else
       IN(stim_2)=0;
   end   
end

%

if t*dt>ts3_start
   if t*dt<ts3_end
       IN(stim_3)=10;
   else
       IN(stim_3)=0;
   end   
end

if t*dt>ts4_start
   if t*dt<ts4_end
       IN(stim_4)=10;
   else
       IN(stim_4)=0;
   end   
end

%}


% Synaptic variable
I_EE=(-I_EE/tau_EE + Df.*sum(S_EE(:,fired),2) )*dt + I_EE;

% Voltage
dVdt=1/taum*(-R.*(V-VL) +IN +I_EE');
V=dVdt*dt + V;


% Representative cell
V_repr(t)=V(repr);
I_repr_ext(t)=IN(repr);
I_repr_EE(t)=I_EE(repr);

if t>1
if V_repr(t)>=VT
    V_repr(t-1)=10;     % add stick to the repr neuron
end
end

% FIRINGS Proccessing
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

%%

%% FINAL PLOT / ONLINE PLOT

if (t-1)*dt==Tframe
    
subplot(2,2,1);
if isempty(firings)==0
plot(firings(:,1)*dt,firings(:,2),'.','MarkerSize',5);
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
imagesc((reshape((IN +I_EE'),sqrt(N),sqrt(N)))',[0 gEE_mean*4]); % Voltage distribution
set(gca,'Ydir','normal');
set(gca,'FontSize',20);             % set the axis with big font
colorbar;
title('Synaptic input (\muA/cm^2)');
box off;

%%

MOV(frame)=getframe(gcf);

frame=frame+1;          % counter for the movie
Tframe=Tframe + dTframe;
%%

end

end

% Save the movie
%movie2avi(MOV,'LIF_net.avi','fps',10,'quality',1);

