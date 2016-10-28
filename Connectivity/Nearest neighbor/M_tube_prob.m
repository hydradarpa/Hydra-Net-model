function A = M_tube_prob(n,m,p)

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
p_connect = makedist('Multinomial','Probabilities',p);
p_rand_all=zeros(1,m*n);        % connections per element

for x=1:1:m         % loop over all all graph edges
    for y=1:1:n

    % number of connections per given neuron
    p_rand=random(p_connect);
    p_rand_all(graph_element(x,y,m))=p_rand;
    
    % connections per given neuron
    connections=zeros(4,1);
    connections(1:p_rand)=1;
    connections=datasample(connections,4,'Replace',false);
            
    % tube border conditions (continous border conditions)
        
        A(graph_element(x,1,m),graph_element(x,n,m))=connections(1);                                    
        A(graph_element(x,n,m),graph_element(x,1,m))=connections(1);
        
        A(graph_element(m,y,m),graph_element(1,y,m))=0;        
        A(graph_element(1,y,m),graph_element(m,y,m))=0;
         
        
    % formula for the rest of the elements
        
       if x>1
       A(graph_element(x,y,m),graph_element(x-1,y,m))=connections(1);            
       A(graph_element(x-1,y,m),graph_element(x,y,m))=connections(1);            
       end       
   
       if y>1      
        A(graph_element(x,y,m),graph_element(x,y-1,m))=connections(2);            
        A(graph_element(x,y-1,m),graph_element(x,y,m))=connections(2);            
       end                                
           
       if x+1<=m
       A(graph_element(x,y,m),graph_element(x+1,y,m))=connections(3);            
       A(graph_element(x+1,y,m),graph_element(x,y,m))=connections(3);            
       end
      
       if y+1<=n
       A(graph_element(x,y,m),graph_element(x,y+1,m))=connections(4);            
       A(graph_element(x,y+1,m),graph_element(x,y,m))=connections(4);            
       end
                            
       
    end
end


% Proof-reading, write more adequate algorithm later
%
B=A;                 % replace the connectivity
        
for i=1:1:m*n         % loop over all all element, half of the symmetric matrix
    
    connected=find(B(i,:)>0);   % indexes of connected elements
    p_rand=p_rand_all(i);   
    
    % remove connections
    if length(connected) > p_rand % cut elements if there too many connections        
    connection_cut=datasample(connected,abs(length(connected)-p_rand),'Replace',false);         
    B(i,connection_cut)=0;
    B(connection_cut,i)=0;
    end
                  
end

A=B;
%}


end

