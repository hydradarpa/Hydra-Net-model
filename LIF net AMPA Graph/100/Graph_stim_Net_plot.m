
%% Generates the plot of the graph with the input stimuli PROBLEM - NOT ALL FIRED CELLS ARE SHOWN!!!

figure('units','normalized','outerposition',[0 0 0.7 0.7]); % show figure window

A_graph=graph(A);

%h=plot(A_graph,'Layout','force')  % force graph

[X,Y]=meshgrid(1:20,1:5);        
x=reshape(X,[1,100]);
y=reshape(Y,[1,100]);
h=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
%xlim([-6 20]);
%ylim([0 53]);

box off

%highlight(h,firings(:,2),'Markersize',7,'Nodecolor','r')    % firing ensemble
%highlight(h,stim_1,'Markersize',7,'Nodecolor','g')          % stimulated cells

title('Green: Stimulated, Red: Fired')
set(gca,'Fontsize',20)
%%


%% Generates the graph with the control nodes
figure('units','normalized','outerposition',[0 0 0.7 0.7]); % show figure window
A_graph=graph(A);


% Plot on tube
%{
[X,Y]=meshgrid(1:13,1:52);        
x=reshape(X,[1,676]);
y=reshape(Y,[1,676]);
h=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
xlim([-6 20]);
ylim([0 53]);
%}

% Force stabilized graph
h=plot(A_graph,'Layout','force')  % force graph

box off

highlight(h,drivernodes,'Markersize',7,'Nodecolor','r')    % control nodes
title('Control nodes')
set(gca,'Fontsize',20)

%%



%%  Play the movie with spikes on the network
% Generates a movie of cells activated during stimulation

firings_TN=zeros(N,t);    % all frings in binary format
                          % neuron ID x time step

for i=1:length(firings)    
    firings_TN(firings(i,2),firings(i,1))=1;
end


figure('units','normalized','outerposition',[0 0 0.7 0.7]); % show figure window

A_graph=graph(A);

%h=plot(A_graph,'Layout','Force')       %handle for the graph plot

%
[X,Y]=meshgrid(1:13,1:52);        
x=reshape(X,[1,676]);
y=reshape(Y,[1,676]);
h=plot(A_graph,'Xdata',x,'Ydata',y);    %handle for the graph plot
%}

hold on

for i=1:10:t
    
    if i*dt>ts1_start
        if i*dt<ts1_end
        highlight(h,stim_1,'Markersize',10,'Nodecolor','g')          % stimulated cells during stimulation  
        else
        highlight(h,stim_1,'Markersize',5,'Nodecolor','g')          % stimulated cells after stimulation
        end
    end
        
    
    
    neuron_fire=find(firings_TN(:,i)>0);
    if isempty(neuron_fire)==0    
    highlight(h,neuron_fire,'Markersize',10,'Nodecolor','r')    % highlight ensemble
    else
    % highlight(h,(1: N),'Markersize',2,'Nodecolor','b')    % highlight ensemble
    end
        
    title(sprintf('Time = %d ms', round(i*dt)));      
    xlim([-6 20]);
    ylim([0 53]);
    
    pause(0.003);
    
   % MOV(t)=getframe(gcf);    
    %savefig(sprintf('Graph_stim_net %d ms .fig', i*dt));    
    
end

%%


%% plot the raster and simulated network structure

figure('units','normalized','outerposition',[0 0 0.7 0.7]); % show figure window
A_graph=graph(A);

subplot(1,2,2)
plot(firings(:,1)*dt,firings(:,2),'.','MarkerSize',10,'Color','b');
hold on;
idx_stim_1 = find(ismember(firings(:,2), burster_idx));         % finds index of stim_1 in firings array
plot(firings(idx_stim_1,1)*dt,firings(idx_stim_1,2),'.','MarkerSize',10,'Color','r');

legend('Excitable','Bursting');


ylabel('Cell index');
xlabel('Time (ms)')
ylim([0 N]);
set(gca,'FontSize',20);             % set the axis with big font
title(sprintf('LIF population, T=%d ms',round(Tframe)));
set(gca,'FontSize',20);             % set the axis with big font
box off;


subplot(1,2,1)

%h=plot(A_graph,'Layout','force')  % forced graph

%
[X,Y]=meshgrid(1:20,1:5);          % plane graph
x=reshape(X,[1,100]);
y=reshape(Y,[1,100]);
h=plot(A_graph,'Xdata',x,'Ydata',y);  %handle for the graph plot
%}

box off
highlight(h,firings(:,2),'Markersize',7,'Nodecolor','b')    % firing ensemble
highlight(h,burster_idx,'Markersize',7,'Nodecolor','r')          % stimulated cells
title('Network structure')
legend('Excitable','Bursting');
%xlim([-6 20]);
%ylim([0 53]);
set(gca,'Fontsize',20)

%%


%% Stimulus

plot((1:1:t)*dt,IN_stim);

xlabel('Time (ms)');
ylabel('Input (pA)');
set(gca,'FontSize',20);             % set the axis with big font
box off;


%%

