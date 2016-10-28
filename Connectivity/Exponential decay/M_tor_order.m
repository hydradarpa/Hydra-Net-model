function A = M_tor_order(n,m,max_order)

% CREATES the adjacency matrix of a size N for elements on the thor (4 neighbours)
% lattice with periodic / empty border conditions

function ind = graph_element(x,y,m)
%UNTITLED2 Summary of this function goes here
%   returns the number of an element on the graph, numbrering: 
% left -> right, up -> down
ind=x+(y-1)*m;
end

% N - size parameter, number of elements in the matrix - sqrt(elements)


A=zeros(m*n,m*n); % size of adjacency matrix

for order=1:max_order

for x=1:1:m         % loop over all all graph edges
    for y=1:1:n
                         
    % tube long border conditions
        A(graph_element(x,order,m),graph_element(x,n,m))=1;        
        A(graph_element(x,n,m),graph_element(x,order,m))=1;                
                               
        
        if x>order                                                     % diagonal border
            A(graph_element(x,order,m),graph_element(x-order,n,m))=1;        
            A(graph_element(x-order,n,m),graph_element(x,order,m))=1;
            
            A(graph_element(x,n,m),graph_element(x-order,order,m))=1;                
            A(graph_element(x-order,n,m),graph_element(x,order,m))=1;
        end
        
        if x+order<=m                                                  % diagonal border
            A(graph_element(x,order,m),graph_element(x+order,n,m))=1;        
            A(graph_element(x+order,n,m),graph_element(x,order,m))=1;                
            
            A(graph_element(x,n,m),graph_element(x+order,order,m))=1;                
            A(graph_element(x+order,n,m),graph_element(x,order,m))=1;            
        end
        
    
   % tube short border conditions    
        A(graph_element(order,y,m),graph_element(m,y,m))=1;                        
        A(graph_element(m,y,m),graph_element(order,y,m))=1;

        
        if y>order                                                     % diagonal border
            A(graph_element(order,y,m),graph_element(m,y-order,m))=1;                        
            A(graph_element(m,y-order,m),graph_element(order,y,m))=1;
            
            A(graph_element(m,y,m),graph_element(order,y-order,m))=1;
            A(graph_element(order,y-order,m),graph_element(m,y,m))=1;
        end
        
        if y+order<=n                                                  % diagonal border
            A(graph_element(order,y,m),graph_element(m,y+order,m))=1;                        
            A(graph_element(m,y+order,m),graph_element(order,y,m))=1;
            
            A(graph_element(m,y,m),graph_element(order,y+order,m))=1;
            A(graph_element(order,y+order,m),graph_element(m,y,m))=1;
        end
        
        
   % maximal diagonal elements, set to zero, physically plausible
   
   A(graph_element(order,order,m),graph_element(m,n,m))=0;
   A(graph_element(m,n,m),graph_element(order,order,m))=0;
   
   A(graph_element(m,order,m),graph_element(order,n,m))=0;
   A(graph_element(order,n,m),graph_element(m,order,m))=0;
  
   
    % formula for the rest of the elements
        
       if x>order
       A(graph_element(x,y,m),graph_element(x-order,y,m))=1;            
       A(graph_element(x-order,y,m),graph_element(x,y,m))=1;
        if y>order                                                      % diagonal elements
            A(graph_element(x,y,m),graph_element(x-order,y-order,m))=1;              
            A(graph_element(x-order,y-order,m),graph_element(x,y,m))=1;
        end
        if y+order<=n                                                   % diagonal elements
            A(graph_element(x,y,m),graph_element(x-order,y+order,m))=1;            
            A(graph_element(x-order,y+order,m),graph_element(x,y,m))=1;
        end
       end       
   
       
       if y>order      
        A(graph_element(x,y,m),graph_element(x,y-order,m))=1;
        A(graph_element(x,y-order,m),graph_element(x,y,m))=1;        
       end                                
           
       
       if x+order<=m
       A(graph_element(x,y,m),graph_element(x+order,y,m))=1;                   
       A(graph_element(x+order,y,m),graph_element(x,y,m))=1;
           if y>order                                                   % diagonal elements
            A(graph_element(x,y,m),graph_element(x+order,y-order,m))=1;            
            A(graph_element(x+order,y-order,m),graph_element(x,y,m))=1;
           end                                                     
            if y+order<=n                                               % diagonal elements
            A(graph_element(x,y,m),graph_element(x+order,y+order,m))=1;            
            A(graph_element(x+order,y+order,m),graph_element(x,y,m))=1;
            end
       end
      
       if y+order<=n
       A(graph_element(x,y,m),graph_element(x,y+order,m))=1;
       A(graph_element(x,y+order,m),graph_element(x,y,m))=1;              
       end
                            
       
    end
end

end

end

