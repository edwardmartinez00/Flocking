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
    
%     % calculate point on boundary in front of bird (point of impact)
%     [X_star,dt_star] = sqboundaryImpact(X,V,P);
%     
%     % Allocating space
%     X_starstar = zeros(P.N,2);
    
    for i=1:P.N     % Bird loop
        % Checking if bounds violated
%         if (X_new(i,1)>P.L) || (X_new(i,1)<-P.L) || ... 
%             (X_new(i,2)>P.L) || (X_new(i,2)<-P.L)
%         
%             if (X_star(i,1)>=X_star(i,2)) && (X_star(i,1)>-X_star(i,2))
%                 normal = [1 0];     % East wall
%             end
%             if (X_star(i,1)<=-X_star(i,2)) && (X_star(i,1)>X_star(i,2))
%                 normal = [0 -1];    % South wall
%             end
%             if (X_star(i,1)<=X_star(i,2)) && (X_star(i,1)<-X_star(i,2))
%                 normal = [-1 0];    % West wall
%             end
%             if (X_star(i,1)>=-X_star(i,2)) && (X_star(i,1)<X_star(i,2))
%                 normal = [0 1];     % North wall
%             end
%             
%             V_new(i,:) = V(i,:)-absorbancy*dot(V(i,:),normal)*normal;
%             X_new(i,:) = X_star(i,:)+(P.dt-dt_star(i,:)).*V_new(i,:);
%             
%             [~,dt_starstar] = sqboundaryImpact(X,V,P);
%             if dt_starstar(i,:) < (P.dt-dt_star)
%                 v_star = V_new(i,:);
%                 x_starstar = X_star(i,:) + (dt_starstar(i,:)-dt_star(i,:)).*V_new(i,:);
%             
%                 if (x_starstar(1)>=x_starstar(2)) && (x_starstar(1)>-x_starstar(2))
%                     normal = [1 0];     % East wall
%                 elseif (x_starstar(1)<=-x_starstar(2)) && (x_starstar(1)>x_starstar(2))
%                     normal = [0 -1];    % South wall
%                 elseif (x_starstar(1)<=x_starstar(2)) && (x_starstar(1)<-x_starstar(2))
%                     normal = [-1 0];    % West wall
%                 else
%                     normal = [0 1];     % North wall
%                 end
%                 V_new(i,:) = v_star-absorbancy*dot(v_star,normal)*normal;
%                 X_new(i,:) = x_starstar + (P.dt-dt_starstar(i,:)).*V_new(i,:);
%             end
%         end  
        if X_new(i,1) > P.L
            dt_star = (P.L-X(i,1))/(V(i,1)+1e-10);
            normal = [1,0];
            change = 1;
        elseif X_new(i,2) > P.L
            dt_star = (P.L-X(i,2))/(V(i,2)+1e-10);
            normal = [0,1];
            change = 1;
        elseif X_new(i,1) < -P.L
            dt_star = (-P.L-X(i,1))/(V(i,1)+1e-10);
            normal = [-1,0];
            change = 1;
        elseif X_new(i,2) < -P.L 
            dt_star = (-P.L-X(i,2))/(V(i,2)+1e-10);
            normal = [0,-1];
            change = 1;
        else
            change = 0;
        end
                    %PUT IN AGAIN
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
