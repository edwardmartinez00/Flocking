%% [X,V] = squareReflex(X,V,X_new,V_new,P)
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
    
% calculate point on boundary in front of bird (point of impact)
[X_star,dt_star] = sqboundaryImpact(X,V,P);
    
    for i=1:P.N     % Bird loop
        % Checking if bounds violated
        if (X_new(i,1)>P.L) || (X_new(i,1)<-P.L) || ... 
            (X_new(i,2)>P.L) || (X_new(i,2)<-P.L)
        
            if (X_star(i,1)>=X_star(i,2)) && (X_star(i,1)>-X_star(i,2))
                normal = [1 0];     % East wall
            elseif (X_star(i,1)<=-X_star(i,2)) && (X_star(i,1)>X_star(i,2))
                normal = [0 -1];    % South wall
            elseif (X_star(i,1)<=X_star(i,2)) && (X_star(i,1)<-X_star(i,2))
                normal = [-1 0];    % West wall
            else
                normal = [0 1];     % North wall
            end
            
            V_new(i,:) = V(i,:)-absorbancy*dot(V(i,:),normal)*normal;
            X_new(i,:) = X_star(i,:)+(P.dt-dt_star(i,:)).*V_new(i,:); 
        end  
%         if X_new(i,1) > P.L
%             dt_star = (P.L-X(i,1))/V(i,1);
%             normal = [1,0];
%             change = 1;
%         elseif X_new(i,2) > P.L
%             dt_star = (P.L-X(i,2))/V(i,2);
%             normal = [0,1];
%             change = 1;
%         elseif X_new(i,1) < -P.L
%             dt_star = (-P.L-X(i,1))/V(i,1);
%             normal = [-1,0];
%             change = 1;
%         elseif X_new(i,2) < -P.L 
%             dt_star = (-P.L-X(i,2))/V(i,2);
%             normal = [0,-1];
%             change = 1;
%         else
%             change = 0;
%         end
%         
%         % Updating position and velocity
%         if change==1 % Outside boundary
%             V_new(i,:) = V(i,:)-absorbancy*dot(V(i,:),normal)*normal;
%             x_star = X(i,:)+dt_star*V(i,:);
%             X_new(i,:) = x_star+(P.dt-dt_star)*V_new(i,:);
%         end
    end
end
