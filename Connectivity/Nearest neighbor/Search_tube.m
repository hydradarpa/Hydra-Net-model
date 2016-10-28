

n=10;
m=6;

% create connectivity matrix
A=M_tube_prob_cut(n,m,p);              

% search for tube with all connected elements, no zero connections

while min(sum(A))==0               
    
    A=M_tube_prob_cut(n,m,p);              % create connectivity matrix
    
end

A_graph=graph(A);
plot(A_graph);

%figure;
%connect_test(A);