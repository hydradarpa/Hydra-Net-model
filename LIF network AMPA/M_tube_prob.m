function A = M_tube_prob(n,m)

% CREATES the adjacency matrix of a size n by m, where n - rows, m -
% columns

% CUTS the connections to fit the probability distrubution p.
% p - vector of the outcomes consisting of probabilites

function ind = graph_element(x,y,m)
%UNTITLED2 Summary of this function goes here
%   returns the number of an element on the graph, numbrering: 
% left -> right, up -> down
ind=x+(y-1)*m;
end

% N - size parameter, number of elements in the matrix - sqrt(elements)

A=zeros(m*n,m*n); % size of adjacency matrix

for x=1:1:m         % loop over all all graph edges
    for y=1:1:n
                         
    % tube border conditions
    % y and y borders
        A(graph_element(x,1,m),graph_element(x,n,m))=1;        
        A(graph_element(1,y,m),graph_element(m,y,m))=0;
    
    % symmetric borders
        A(graph_element(x,n,m),graph_element(x,1,m))=1;        
        A(graph_element(m,y,m),graph_element(1,y,m))=0;        
        
       if x>1
       A(graph_element(x,y,m),graph_element(x-1,y,m))=1;            
       A(graph_element(x-1,y,m),graph_element(x,y,m))=1;              
       end       
    
       if y>1      
        A(graph_element(x,y,m),graph_element(x,y-1,m))=1;
        A(graph_element(x,y-1,m),graph_element(x,y,m))=1; 
       end                                
    
       % formula for the rest of the elements
       if x+1<=m
       A(graph_element(x,y,m),graph_element(x+1,y,m))=1;                   
       A(graph_element(x+1,y,m),graph_element(x,y,m))=1;                     
       end
      
       if y+1<=n
       A(graph_element(x,y,m),graph_element(x,y+1,m))=1;
       A(graph_element(x,y+1,m),graph_element(x,y,m))=1;              
       end
                            
    end
end


% Cut off connections with given probabilities

%p_conect = makedist('Multinomial','Probabilities',p);

%for i=1:1:m*n         % loop over all all elements
    
%    connected=find(A(i,:)>0); % indexes of connected elements
%    if length(connected) ~= random(p_conect)    % compare with dist
%        connection_cut=randperm(connected,abs(connected-random(p_connect)));
    
%end

end

