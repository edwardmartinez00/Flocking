%% V_new = squareImpact(X,V,P)
%
% Checks points at a soft square boundary using only impact as a reference
%
% X, position vector
% V, velocity vector
% P.N, number of birds
% P.L, apothem of square
% P.d, distance that W occurs in
% P.dt, time step size
%
% BoundaryX, vector of 4 conditions for boundary of X
% BoundaryY, vector of 4 conditions corresponding to boundary of Y
function [X,V,BoundaryX,BoundaryY] = squareSoft(X,V,P)
    V_functions;

    for i = 1:P.N
        % Time to cross right(left...) boundary
        dtR = (P.L-X(i,1))/V(i,1);
        dtL = (-P.L-X(i,1))/V(i,1);
        dtU = (P.L-X(i,2))/V(i,2);
        dtD = (-P.L-X(i,2))/V(i,2);
        impactTime = [dtR dtL dtU dtD];
    
        % Finding boundary affecting vectors
        M = max(impactTime);
        impactTime(impactTime < 0) = M+1; % Set negative vectors large
        dt_star = min(impactTime);      % Time till impact closest wall
    
        % Calculating new X and V
        x_star = X(i,:) + dt_star*V(i,:);   % Wouldbe impact point
        diffstar = x_star-X(i,:);
        r = norm(diffstar);     % Distance of bird to impact
        v_new = V(i,:) - P.dt*w(r,P)*(diffstar/r);  % Euler's Method
        X(i,:) = x_star + (P.dt-dt_star)*V(i,:);
        V(i,:) = v_new;     % Update V
    end
BoundaryX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
BoundaryY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)];             
end