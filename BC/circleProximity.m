%% dV = circleProximity(X,V,P)
%
% Finds change in velocity at a soft circle boundary based solely on the
% closest point to the birds
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position depending on distance to boundary
% P.N, number of birds
% P.R, radius of circle
% P.dt, time step size
function dV = circleProximity(X,V,P)
    V_functions;

    dV = noBounds(X,V,P);
    
    % allocating space
    X_hat = zeros(P.N,2);
    
    for i = 1:P.N     % Bird loop
        X_hat(i,:) = P.R*X(i,:)/norm(X(i,:));        % closest wall
    end   
    
    % Calculating distance to closest wall
    D_hat = X_hat-X;
    r_hat = sqrt(sum(D_hat.^2,2));
    r_hat = r_hat*[1 1];
    
    dV = dV - w(r_hat,P).*(D_hat./r_hat);
end