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
function [X,V,BoundaryX,BoundaryY] = circleImpact(X,V,X_new,P)
    V_functions;

    for i = 1:P.N     % Bird loop
        if norm(X_new(i,:)) > P.R-P.d % Boundary effect
            a = norm(V(i,:))^2;
            b = 2*dot(X(i,:),V(i,:));
            c = (norm(X(i,:)))^2-P.R^2;
            dt_star = (-b+sqrt(b^2-4*a*c))/(2*a);

            x_star = X(i,:) + dt_star*V(i,:);
            x_hat = P.R*X(i,:)/norm(X(i,:));

            hatdiff = x_hat-X(i,:);
            r = norm(X(i,:)-x_star);

            X(i,:) = x_star + (P.dt-dt_star)*V(i,:);
            V(i,:) = V(i,:) - P.dt*w(r,P)*(hatdiff)/norm(hatdiff);

        else %if does not feel boundary
            X(i,:) = X(i,:) + P.dt*V(i,:);
        end
    end
theta=linspace(0,2*pi);
BoundaryX=P.R*cos(theta);
BoundaryY=P.R*sin(theta);  
end