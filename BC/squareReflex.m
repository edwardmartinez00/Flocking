%% [X,V,BoundaryX,BoundaryY] = squareReflex(X,V,X_new,dV,P)
%
% Checks points at a reflexive square boundary
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position vector (depending on relation to bounds)
% P.N, number of birds
% P.L, apothem of square
% P.absorb, boundary absorb proportion
% P.dt, time step size
% dV, change in V after flocking functions
%
% BoundaryX, vector of 4 conditions for boundary of X
% BoundaryY, vector of 4 conditions corresponding to boundary of Y
function [X,V,BoundaryX,BoundaryY] = squareReflex(X,V,X_new,dV,P)
    absorbancy = 2 - P.absorb;
    
    for i=1:P.N     % Bird loop
        % Checking if bounds violated
        if X_new(i,1) > P.L
            dt_star = (P.L-X(i,1))/V(i,1);
            normal = [1,0];
            change = 1;
        elseif X_new(i,2) > P.L
            dt_star = (P.L-X(i,2))/V(i,2);
            normal = [0,1];
            change = 1;
        elseif X_new(i,1) < -P.L
            dt_star = (-P.L-X(i,1))/V(i,1);
            normal = [-1,0];
            change = 1;
        elseif X_new(i,2) < -P.L 
            dt_star = (-P.L-X(i,2))/V(i,2);
            normal = [0,-1];
            change = 1;
        else
            change = 0;
        end
        
        % Updating position and velocity
        if change==1 % Outside boundary
            v_new = V(i,:)-absorbancy*dot(V(i,:),normal)*normal;
            X_star = X(i,:)+dt_star*V(i,:);
            X(i,:) = X_star+(P.dt-dt_star)*v_new;
            V(i,:) = v_new;
        else % Inside boundary
            X(i,:)=X(i,:)+P.dt*V(i,:);
        end
    end
    V = V + P.dt*dV;    % update V with flocking
end
