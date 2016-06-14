%% dV = squareImpact(X,V,P)
%
% Finds change in velocity at a soft square boundary solely based on impact 
% as a reference
%
% X, position vector
% V, velocity vector
% P.N, number of birds
% P.L, apothem of square
% P.d, distance that W occurs in
% P.dt, time step size
function dV = squareImpact(X,V,P)
    V_functions;

    % calculate acceleration without boundary conditions
    dV = noBounds(X,V,P);
    
    % point on boundary in front of birds
    X_star = sqboundaryImpact(X,V,P);        % would be impact point
    
    % Calculating distance to impact
    D_star = X_star-X;
    r_star = sqrt(sum(D_star.^2,2));
    r_star = r_star*[1 1];    
    
    dV = dV - w(r_star,P).*(D_star./r_star); 
end