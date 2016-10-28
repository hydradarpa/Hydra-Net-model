
%% Generate the graph

%A=M_tube_prob_cut_add(40,40,p);    % 4 potential connections

n=13;
m=52;

A=M_FullTube_prob_cut_add(n,m,p(1:4),3); % 8 potential connections
A_graph=graph(A);

%% Plot the graph
figure

subplot(1,2,1)
plot(A_graph,'Layout','force');
title(sprintf('Graph plot (%d x %d nodes)',n,m));
box off;
set(gca,'Fontsize',20);

subplot(1,2,2)
histogram(sum(A));
ylabel('Number of neurons');
xlabel('Number of connections');
title('Max = 4 connections');
box off;
set(gca,'Fontsize',20);

%% Adjacency matrix visualization

figure

imagesc(A);
xlabel('Neuron ID');
ylabel('Neuron ID');
title(sprintf('Adjacency matrix (%d x %d nodes)',n,m));

set(gca,'Fontsize',20);

%%



%% Plot the connection probability and compare with data
figure

A_sum=sum(A);              % take summ of all the elements

A_conn=[length(find(A_sum==1)),length(find(A_sum==2)),...
   length(find(A_sum==3)),length(find(A_sum==4)),...
   length(find(A_sum==5)), length(find(A_sum==6)),...
   length(find(A_sum==7)), length(find(A_sum==8)),...
   length(find(A_sum==9)),length(find(A_sum==10))];

% Connectivity distribution
A_conn=A_conn/length(A);
A_conn=A_conn(find(A_conn>0));

subplot(1,2,1)
bar(A_conn);
ylabel('Connection probability');
xlabel('Number of connections');
ylim([0 1])
title('Model')
box off;
set(gca,'Fontsize',20);

subplot(1,2,2)
bar(p(find(p>0)));
ylabel('Connection probability');
xlabel('Number of connections');
title('Data');
ylim([0 1])
box off;
set(gca,'Fontsize',20);

%%