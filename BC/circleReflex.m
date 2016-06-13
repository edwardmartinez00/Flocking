%% [X,V,BoundaryX,BoundaryY] = circleReflex(X,V,X_new,P)
%
% Checks points at a reflexive/absorbant circle boundary
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position vector (depending on relation to bounds)
% P.N, number of birds
% P.R, radius of circle
% P.absorb, boundary absorb proportion
% P.dt, time step size
%
% BoundaryX, vector of 4 conditions for boundary of X
% BoundaryY, vector of 4 conditions corresponding to boundary of Y
function [X,V] = circleReflex(X,V,X_new,dV,P)
    absorbancy = 2 - P.absorb;
    
    for i = 1:P.N     % Bird loop
        % Dirty Code
        if norm(X(i,:)) > P.R 
            X(i,:) = P.R*X(i,:)/norm(X(i,:));
            X_new = X + P.dt*V;
        end
        
        if norm(X_new(i,:)) > P.R    % Outside bound
            a = norm(V(i,:))^2;
            b = 2*dot(X(i,:),V(i,:));
            c = (norm(X(i,:)))^2-P.R^2;
            
            dt_star = (-b + sqrt(b^2-4*a*c))/(2*a);     % time of impact
            x_star = X(i,:) + dt_star*V(i,:);           % point of impact     
            normal = X(i,:)/norm(X(i,:));               % normal vector
            
            V(i,:) = V(i,:) - absorbancy*dot(V(i,:),normal)*normal;
            X(i,:) = x_star + (P.dt-dt_star)*V(i,:);
            
%             if X(i,:) > P.R         % if still outside, put on boundary
%                 X(i,:) = P.R*X(i,:)/norm(X(i,:));
%             end
        else
%             if X(i,:) > P.R         % if still outside, put on boundary
%                 X(i,:) = P.R*X(i,:)/norm(X(i,:));
%             end
            X(i,:) = X_new(i,:);
        end
    end
    V = V + P.dt*dV;        % updates V with flocking
end