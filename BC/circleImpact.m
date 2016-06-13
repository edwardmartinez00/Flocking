%% [X,V,BoundaryX,BoundaryY] = circleSoft(X,V,X_new,P)
%
% Checks points at a soft cirlce boundary
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position depending on distance to boundary
% P.N, number of birds
% P.R, radius of circle
% P.dt, time step size
%
% BoundaryX, vector of 4 conditions for boundary of X
% BoundaryY, vector of 4 conditions corresponding to boundary of Y
function dV = circleImpact(X,V,P)
    V_functions;

    dV = noBounds(X,V,P);
    
    % allocating space
    X_star = zeros(P.N,2);
    X_hat = zeros(P.N,2);
    
    for i = 1:P.N     % Bird loop
        a = norm(V(i,:))^2;
        b = 2*dot(X(i,:),V(i,:));
        c = (norm(X(i,:)))^2-P.R^2;

        dt_star = (-b+sqrt(b^2-4*a*c))/(2*a);   % impact time
        X_star(i,:) = X(i,:) + dt_star*V(i,:);       % impact point
        X_hat(i,:) = P.R*X(i,:)/norm(X(i,:));        % closest wall

%         hatdiff = X_hat-X(i,:);
%         r = norm(X(i,:)-X_star);            % distance to impact point

%         X(i,:) = x_star + (P.dt-dt_star)*V(i,:);
%         V(i,:) = V(i,:) - P.dt*w(r,P)*(hatdiff)/norm(hatdiff);
%         dV(i,:) = flock_dV(i,:) - w(r,P)*(hatdiff)/norm(hatdiff);
    end  
    
    % Calculating distance to impact
    D = X_star-X;
    r = sqrt(sum(D.^2,2));
    r = r*[1 1];    
    
    % Calculating distance to closest wall
    D_hat = X_hat-X;
    r_hat = sqrt(sum(D.^2,2));
    r_hat = r_hat*[1 1];
    
    dV = dV - w(r,P).*(D_hat./r_hat);
end