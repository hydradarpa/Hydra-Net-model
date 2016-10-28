function A = exponential_horn(N_x,N_y,p0,sigma_x)
% This function generates the adjacency matrix of the graph N_x by N_y size
% Connection probabilities on this graph are given on the rectangular grid
% by the connection probability formula

%  P(x,y,x0,y0) = p0 x exp[-(x-x0)^2/sigma_x -(y-y0)^2/sigma_y]

% where x, y - location of the element on the graph
% and x0, y0 - location of a given element

% autapses are not allowed (are there any?)


A=zeros(N_x*N_y,N_y*N_x);

% loop over all elements
for i=1:N_x*N_y
    for j=1:N_y*N_x
        
% loop oveer all connections of the element        
%        for x=1:N_x*N_y
%            for y=1:N_y*N_x
                p_random=rand;
                p=p0*exp(-(i-j)^2/sigma_x );    % probability of being connected
                
                %A(i,j)=p;
                
                %  Need to implement the periodic BORDER FOR PROBABILITY
                
                % This part works!
                %
                if p_random<=p
                    A(i,j)=1;      % connections are symmetric CHECK THE CONFLICTS!
                    A(j,i)=1;   
                elseif p_random>p
                    A(i,j)=0;
                    A(j,i)=0;
                end
                %}
                
%            end
%        end
        
% autapses
A(i,i)=0;        

        
    end
    
end


end

