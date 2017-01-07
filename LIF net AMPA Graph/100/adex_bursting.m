
clear;
%close all;

% Integration time and time step
N=500;
dt=0.05;

% cortex-matched model, Brette-Touboul 2008
c=281; 
gl=30;
el=-70.6; 
vt=-50.4;
delta=2; 
vreset=-47.2;
a=4; tauw=500; b=80;      % time constant to determine the rhythm
% spike mark and number of spikes
vspike=0;
nsp=0;

% model variables
v=zeros(1,round(N/dt));
w=zeros(1,round(N/dt));

% stimuli
stimulus=0; % pA
input=stimulus*ones(1,round(N/dt));

% initial conditions
v(1)=-44.15; %Vrest
w(1)=102.83;


for i=2:1:round(N/dt)
        t(i)=(i-1)*dt;
                                                        
        v(i)=dt/c*(-gl*(v(i-1)-el)+gl*delta*exp((v(i-1)-vt)/delta)-w(i-1)+input(i)) + v(i-1);
        w(i)=dt/tauw*(a*(v(i-1)-el)-w(i-1)) + w(i-1);
                                       
        if  v(i)>=vspike
            v(i-1)=0;
            v(i)=vreset;
            w(i)=w(i) + b;
            nsp=nsp + 1;
        end                                                            
                
end 

figure

subplot(2,1,1)
plot(t,v)
box off
set(gca,'Fontsize',20)
ylabel('Voltage (mV)')
title(sprintf('AdEx neuron (%d pA)',input(1)))

subplot(2,1,2)
plot(t,w)
box off
ylabel('Adaptation (pA)')
xlabel('time (ms)')
set(gca,'Fontsize',20)
