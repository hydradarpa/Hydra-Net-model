function [mean_connect, std_connect] = connect_test(A)

% This function takes the adjacency matrix and plots the number of
% connections in each element

close all;

figure;
connections=sum(A);                         % sum of all connections in the matrx
histogram(connections);
set(gca,'Fontsize',30);
ylabel('Count');
xlabel('Number of connections');
title(sprintf('Connectivity test %d nodes',length(A)));
box off

% save the mean and std of the distribution
mean_connect=mean(connections);             % average number of connections
std_connect=std(connections)/length(A);     % standard error

% plotting part
%savefig(sprintf('Connect_test %d nodes.fig',length(A)));           % safe *.fig
%saveas(gcf,sprintf('Connectivity test %d nodes.jpg',length(A)),'jpg');  % save *.jpg


end

