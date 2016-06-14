%% dV = squareSoft(X,V,P)
%
% Finds change in velocity at a soft square boundary with both impact 
% and closest boundary
%
% X, position vector
% V, velocity vector
% P.N, number of birds
% P.L, apothem of square
% P.d, distance that W occurs in
% P.dt, time step size
function dV = squareSoft(X,V,P)
    V_functions;

    % calculate acceleration without boundary conditions
    dV = noBounds(X,V,P);
    
    % point on boundary in front of birds
    X_star = sqboundaryImpact(X,V,P);        % would be impact point
    
    % point on boundary closest to birds
    X_hat = sqboundaryDistance(X,P);    % closest point 
    
    % Calculating distance to impact
    D_star = X_star-X;
    r_star = sqrt(sum(D_star.^2,2));
    r_star = r_star*[1 1];    
    
    % Calculating distance to closest wall
    D_hat = X_hat-X;
    r_hat = sqrt(sum(D_hat.^2,2));
    r_hat = r_hat*[1 1];
    
    dV = dV - w(r_star,P).*(D_hat./r_hat);             
end