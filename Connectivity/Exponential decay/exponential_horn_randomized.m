function A = exponential_horn(N_x,N_y,p0,sigma_x,sigma_y,trials)
% This function generates the adjacency matrix of the graph N_x by N_y size
% Connection probabilities on this graph are given on the rectangular grid
% by the connection probability formula

%  P(x,y,x0,y0) = p0 x exp[-(x-x0)^2/sigma_x -(y-y0)^2/sigma_y]

% where x, y - location of the element on the graph
% and x0, y0 - location of a given element

% autapses are not allowed (are there any?)


% DOES NOT WORK WELL!

A=zeros(N_x*N_y,N_y*N_x);

for tr=1:trials

x_rand=randperm(N_x*N_y);
y_rand=randperm(N_x*N_y);
x_con_rand=randperm(N_x*N_y);
y_con_rand=randperm(N_x*N_y);

% loop over all elements
for i=1:N_x*N_y
    for j=1:N_y*N_x
        
% loop oveer all connections of the element        
        for x=1:N_x*N_y
            for y=1:N_y*N_x
                p_random=rand;
                p=p0*exp(-(x-x_con_rand(i))^2/sigma_x -(y-y_con_rand(j))^2/sigma_y);    % probability of being connected
                
                if p_random<=p
                    A(x,y)=1;      % connections are symmetric CHECK THE CONFLICTS!
                    A(y,x)=1;   
                elseif p_random>p
                    A(x,y)=0;
                    A(y,x)=0;
                end
                
            end
        end
        
% no reccurent synapses
%A(i,j)=0;        

        
    end
    
end


end

end

