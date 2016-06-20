% [X,V] = squareReflex(X,V,X_new,V_new,P)
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
function [X_new,V_new] = squareReflex(X,V,X_new,V_new,P)
    absorbancy = 2 - P.absorb;
    
    for i=1:P.N     % Bird loop
        %Checking if bounds are violated
        if X_new(i,1) > P.L
            dt_star = (P.L-X(i,1))/(V(i,1)+realmin);
            normal = [1,0];
            change = 1;
        elseif X_new(i,2) > P.L
            dt_star = (P.L-X(i,2))/(V(i,2)+realmin);
            normal = [0,1];
            change = 1;
        elseif X_new(i,1) < -P.L
            dt_star = (-P.L-X(i,1))/(V(i,1)+realmin);
            normal = [-1,0];
            change = 1;
        elseif X_new(i,2) < -P.L 
            dt_star = (-P.L-X(i,2))/(V(i,2)+realmin);
            normal = [0,-1];
            change = 1;
        else
            change = 0;
        end
        % Updating position and velocity
        if change==1 % Outside boundary
            V_new(i,:) = V(i,:)-absorbancy*dot(V(i,:),normal)*normal;
            x_star = X(i,:)+dt_star*V(i,:);
            X_new(i,:) = x_star+(P.dt-dt_star)*V_new(i,:);
            
            rounding = 4;
            if round(X_new(i,1),rounding) > P.L
                dt_starstar = (P.L-X(i,1))/(V(i,1)+realmin);
                normal = [1,0];
                change = 1;
            elseif round(X_new(i,2),rounding) > P.L
                dt_starstar = (P.L-X(i,2))/(V(i,2)+realmin);
                normal = [0,1];
                change = 1;
            elseif round(X_new(i,1),rounding) < -P.L
                dt_starstar = (-P.L-X(i,1))/(V(i,1)+realmin);
                normal = [-1,0];
                change = 1;
            elseif round(X_new(i,2),rounding) < -P.L 
                dt_starstar = (-P.L-X(i,2))/(V(i,2)+realmin);
                normal = [0,-1];
                change = 1;
            else
                change = 0;
            end
            
            if change == 1
                v_star = V_new(i,:);
                x_starstar = x_star + (dt_starstar-dt_star).*V_new(i,:);
                V_new(i,:) = v_star-absorbancy*dot(v_star,normal)*normal;
                X_new(i,:) = x_starstar + (P.dt-dt_starstar).*V_new(i,:);
            end
        end
    end
end
