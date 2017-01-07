function A = M_FullTube_smartcut(n,m,p,trials)

% CREATES the adjacency matrix of a size n by m, where n - rows, m -
% columns

% CUTS the connections to fit the probability distrubution p.
% ADD the connections to fit the probability distribution
% p - vector of the outcomes consisting of probabilites


function ind = graph_element(x,y,m)
%UNTITLED2 Summary of this function goes here
%   returns the number of an element on the graph, numbrering: 
% left -> right, up -> down
ind=x+(y-1)*m;
end

function A_conn=A_p_conn(A)        
% Generates the connection probability for the matrix.

A_sum=sum(A);              % take summ of all the elements

A_conn=[length(find(A_sum==1)),length(find(A_sum==2)),...
   length(find(A_sum==3)),length(find(A_sum==4)),...
   length(find(A_sum==5)), length(find(A_sum==6)),...
   length(find(A_sum==7)), length(find(A_sum==8)),...
   length(find(A_sum==9)),length(find(A_sum==10))];

% Connectivity distribution
A_conn=A_conn/length(A);

    
end

% N - size parameter, number of elements in the matrix - sqrt(elements)


A=zeros(m*n,m*n); % size of adjacency matrix


for x=1:1:m         % loop over all all graph edges
    for y=1:1:n
                         
    % tube border conditions    
        A(graph_element(x,1,m),graph_element(x,n,m))=1;        
        A(graph_element(x,n,m),graph_element(x,1,m))=1;                
                               
        
        if x>1                                                     % diagonal border
            A(graph_element(x,1,m),graph_element(x-1,n,m))=1;        
            A(graph_element(x-1,n,m),graph_element(x,1,m))=1;                
        end
        
        if x+1<=m                                                  % diagonal border
            A(graph_element(x,1,m),graph_element(x+1,n,m))=1;        
            A(graph_element(x+1,n,m),graph_element(x,1,m))=1;                
        end
        
        A(graph_element(1,y,m),graph_element(m,y,m))=0;                        
        A(graph_element(m,y,m),graph_element(1,y,m))=0;        
                    
    % formula for the rest of the elements
        
       if x>1
       A(graph_element(x,y,m),graph_element(x-1,y,m))=1;            
       A(graph_element(x-1,y,m),graph_element(x,y,m))=1;
        if y>1                                                      % diagonal elements
            A(graph_element(x,y,m),graph_element(x-1,y-1,m))=1;              
            A(graph_element(x-1,y-1,m),graph_element(x,y,m))=1;
        end
        if y+1<=n                                                   % diagonal elements
            A(graph_element(x,y,m),graph_element(x-1,y+1,m))=1;            
            A(graph_element(x-1,y+1,m),graph_element(x,y,m))=1;
        end
       end       
   
       
       if y>1      
        A(graph_element(x,y,m),graph_element(x,y-1,m))=1;
        A(graph_element(x,y-1,m),graph_element(x,y,m))=1;        
       end                                
           
       
       if x+1<=m
       A(graph_element(x,y,m),graph_element(x+1,y,m))=1;                   
       A(graph_element(x+1,y,m),graph_element(x,y,m))=1;
           if y>1                                                   % diagonal elements
            A(graph_element(x,y,m),graph_element(x+1,y-1,m))=1;            
            A(graph_element(x+1,y-1,m),graph_element(x,y,m))=1;
           end                                                     
            if y+1<=n                                               % diagonal elements
            A(graph_element(x,y,m),graph_element(x+1,y+1,m))=1;            
            A(graph_element(x+1,y+1,m),graph_element(x,y,m))=1;
            end
       end
      
       if y+1<=n
       A(graph_element(x,y,m),graph_element(x,y+1,m))=1;
       A(graph_element(x,y+1,m),graph_element(x,y,m))=1;              
       end
                            
       
    end
end


% Smart cut of connections, should be determined up to 10 elements
p_connect = makedist('Multinomial','Probabilities',p);

B=A;  % matrix we change

for q=1:trials   % run optimization multiple times

horz=randperm(n*m);
vert=randperm(n*m);
    
% Go over all elements in the matrix
for i=1:1:m*n                         % MAKE RANDOM SAMPLING OVER THE ELEMENTS!
   % i=horz(q);
    for j=1:1:m*n
   % j=vert(s);    
                
        if B(i,j)==0            
        else
           B_old=B;                   % matrix on the previos step
           B(i,j)=0;                  % cut the connection
           B(j,i)=0;
           B_graph=graph(B);          % create a graph out of matrix B  
           if max(conncomp(B_graph))==1 % check if all elements are connected               
           
              B_prob=A_p_conn(B);          % connection prob of new matrix
              B_old_prob=A_p_conn(B_old);  % connection prob of old matrix
              
                if sum(abs(B_prob - p)) < sum(abs(B_old_prob - p))   % abs minimization
                   B(i,j)=0;             % cut only if it is beneficial
                   B(j,i)=0;
                else
                    B(i,j)=1;
                    B(j,i)=1;
                end
           else
               B(i,j)=1;
               B(j,i)=1;
           end
           
        end
    end
end

end

A=B;   %

end