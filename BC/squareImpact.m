%% V_new = squareSoft(X,V,P)
%
% Checks points at a soft square boundary
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
function [X,V,BoundaryX,BoundaryY] = squareImpact(X,V,P)
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
    
        % Closest wall
        % Distance from right(left...) boundary
        dR = norm(P.L-X(i,1));
        dL = norm(-P.L-X(i,1));
        dU = norm(P.L-X(i,2));
        dD = norm(-P.L-X(i,2));
        dist = [dR dL dU dD];
        
        % Finding closest wall to bird
        m = min(dist);              % Find closest wall
        dist = dist==m;
        x_hat = dist*[P.L, X(i,2); 
                     -P.L, X(i,2); 
                      X(i,1), P.L;  
                      X(i,1),-P.L];      % Closest point on wall
        diffhat = x_hat - X(i,:);       % Diff of point and closest point on wall
        
        % Calculating new X and V
        x_star = X(i,:) + dt_star*V(i,:);   % Would be impact point
        r = norm(X(i,:)-x_star);     % Distance of bird to impact
        v_new = V(i,:) - P.dt*w(r,P)*(diffhat/norm(diffhat));  % Euler's Method
        X(i,:) = x_star + (P.dt-dt_star)*V(i,:);
        V(i,:) = v_new;     % Update V
    end
BoundaryX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
BoundaryY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)];             
end