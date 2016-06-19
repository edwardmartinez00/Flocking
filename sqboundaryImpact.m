function [X_star,dt_star] = sqboundaryImpact(X,V,P)
% used in the flocking model with soft square potential boundary conditions
% X_star is an Nx2 matrix of the points on the boundary in front of birds
% X is an Nx2 matrix of the current positions of N birds
% V is an Nx2 matrix of the current velocities of N birds
% the domain is the square [-L,L]x[-L,L]

dt_star = zeros(P.N,2);

for i = 1:P.N
    % Time to cross east(west...) boundary
    dt_E = (P.L-X(i,1))/V(i,1);
    dt_W = (-P.L-X(i,1))/V(i,1);
    dt_S = (-P.L-X(i,2))/V(i,2);
    dt_N = (P.L-X(i,2))/V(i,2);   
    impactTime = [dt_E dt_W dt_S dt_N];     % possible impact times
    
    % Finding boundary affecting vectors
    M = max(impactTime);
    impactTime(impactTime < 0) = M+1; % Set negative vectors large
    dt_star(i,:) = min(impactTime);      % Time till impact point in front
end

%t_star = t_star*[1 1];

X_star = X+dt_star.*V;      % Would be impact point
end

