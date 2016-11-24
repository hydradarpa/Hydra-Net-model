%% INTEGRATION PARAMETERS
clear;       % remove previous varibles
%close all;

T=5000;        % total time, mS
dt=0.1;     % time step, ms

%time=(1:1:(T/dt)).*dt;      % TIME VECTOR
Df=1/dt;            % Delta function approximation

Tframe=0;    % movie
dTframe=100; % ms
frame=100;

figure('units','normalized','outerposition',[0 0 0.8 0.8]); % show figure window
%%

%% NETWORK PARAMETERS
N=100;        % Number of neurons
%%

%% CONNECTIVITY
gEE_mean=5;  % mA/cm^2 (current synapses)  # 33.35
gEE_sigma=5;

tau_EE=2;      % ms, AMPA 5.4

load('Adjacency_20X5.mat');


S_EE=A_tube*gEE_mean;    % connectivity matrix

% GAP junctions, the same location as synapses, PROBLEM WITH IMPLEMENTATION
g_GAP=0;      % 0.0001

l=zeros(N,4);

for i=1:1:N 
   l(i,1:length(find(S_EE(i,1:end)>0)))=find(S_EE(i,1:end)>0);
end

% number of connections per element
m=l>0;
m=sum(m,2);

% Representative cell
repr=round(N/2);             % number of representative neuron
%%

%% NEURON PARAMETERS
% cortical neurons, Brette-Gerstner 2005
c(1:N)=281; 
gl(1:N)=random('Normal',30,5,1,N); % heterogenous population
el(1:N)=-70.6; 
vt(1:N)=-50.4;
delta(1:N)=2; 
vreset(1:N)=-70.6;
a(1:N)=4;
tauw(1:N)=144;
b(1:N)=80;

% bursting subpopulation, oscillators are heterogenous
burster_idx=randperm(100,10);
vreset(burster_idx)=-47.2; % random cells
tauw(burster_idx)=random('Normal',500,50,1,length(burster_idx));
%}

% spike mark and number of spikes
vspike=10;
%%

%% STIMULATION
IN(1:N)=0;              % mA/cm^2 to all cells
IN(burster_idx)=800;    % mA/cm^2

% stimulated neuron
stim_1=50;                % central neuron

ts1_start=0;        % start of first stimuli
ts1_end=0;          % end of first stimuli

%%

%% ICs
%Constant ICs
V(1:N)=el;   % mV
W(1:N)=0;
V_repr=zeros(1); % mV
W_repr=zeros(1); % pA

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
%{
if t*dt>ts1_start
    if t*dt<ts1_end
       IN(stim_1)=800;
    else
        IN(stim_1)=0;
    end
end
%}

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

% GAP junciton currents
I_GAP=g_GAP*(m.*V'-sum(repmat(V',1,N).*S_EE,1)');


% Voltage
V=(dt./c).*(-gl.*(V-el)+gl.*delta.*exp((V-vt)./delta)-W +IN +I_EE'+I_GAP') + V;

% FIRINGS Proccessing, avoid W jumps!
fired=[];
fired=find(V>=vspike);

% reset condition       
if isempty(fired)==0
    V(fired)=vreset(fired);
    W(fired)=W(fired)+b(fired);
    else
    W=(dt./tauw).*(a.*(V-el)-W) + W;
    end

% Representative cell
V_repr(t)=V(repr);
W_repr(t)=W(repr);
I_repr_ext(t)=IN(repr);
I_repr_EE(t)=I_EE(repr);
% add stick to the repr neuron
if V_repr(t)>=-40
    V_repr(t-1)=10;
    V_repr(t)=vreset(repr);
end


% Record the spikes
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
plot((1:t)*dt,V_repr,(1:t)*dt,W_repr);
ylabel('mV and pA');
set(gca,'FontSize',20);             % set the axis with big font
title('Representative neuron');
set(gca,'FontSize',20);             % set the axis with big font
legend('Voltage (mV)', 'W (pA)')
box off;

subplot(2,2,3);
imagesc((reshape(V,sqrt(N),sqrt(N)))',[-75 0]); % Voltage distribution
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

