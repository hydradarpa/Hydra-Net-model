function [all_nodes, h_freq] = connect_test_multiple(m,n,p,trials,bins)

% This function takes the adjacency matrix and plots the number of
% connections in each element and loops over realisations in trials

close all;

all_nodes=zeros(trials, m*n);
%h_freq=zeros(trials,bins);
A=zeros(m*n,m*n);

for i=1:1:trials
   
    A=M_tube_prob(m,n,p);              % create connectivity matrix
    
    while min(sum(A))==0               % make sure we take the right matrix
    A=M_tube_prob(m,n,p);              % create connectivity matrix
    end
   
   connections=sum(A);                  % summ all connections in the matrix
   h=histogram(connections,bins);       % create a histogram
%   h_freq(i,:)=h.Values;                % matrix of frequencies
   all_nodes(i,:)=connections;          % save all particular connections
end

%
%errorbar(mean(h_freq),std(h_freq),'Linewidth',2);
histogram(A_all);
set(gca,'Fontsize',30);
ylabel('Number of nodes');
xlabel('Number of connections');
title(sprintf('Connectivity: %d nodes',length(A)));
box off

%}

% plotting part
%savefig(sprintf('Connect_test %d nodes.fig',length(A)));           % safe *.fig
%saveas(gcf,sprintf('Connectivity test %d nodes.jpg',length(A)),'jpg');  % save *.jpg

end

