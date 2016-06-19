%% dV = circleReflex(X,V,P)
%
% Finds change in velocity at a reflexive/absorbant circle boundary
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position vector (depending on relation to bounds)
% P.N, number of birds
% P.R, radius of circle
% P.absorb, boundary absorb proportion
% P.dt, time step size
function [X_new,V_new] = circleReflex(X,V,X_new,V_new,P)
    absorbancy = 2 - P.absorb;
    
    for i = 1:P.N     % Bird loop
        if norm(X_new(i,:)) > P.R    % Outside bound
            a = norm(V(i,:))^2;
            b = 2*dot(X(i,:),V(i,:));
            c = (norm(X(i,:)))^2-P.R^2;
            
            dt_star = (-b + sqrt(b^2-4*a*c))/(2*a);     % time of impact
            x_star = X(i,:) + dt_star*V(i,:);           % point of impact     
            normal = X(i,:)/norm(X(i,:));               % normal vector
            
            V_new(i,:) = V(i,:) - absorbancy*dot(V(i,:),normal)*normal; % change only the things in new
            X_new(i,:) = x_star + (P.dt-dt_star)*V_new(i,:);
        end
    end
    
    %%----------------For Sending Warnings------------%%
%     for i = 1:P.N
%         if norm(X(i,:))>P.R
%             warning('Bird %i Outside',i);
%         end
%     end
end