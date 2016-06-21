function [X,V] = circlePortal(X,V,X_new,V_new,P)
% teleports birds
% Detailed explanation goes here

for i = 1:P.N     % Bird loop
    if norm(X_new(i,:)) > P.R    % Outside bound
        a = norm(V(i,:))^2;
        b = 2*dot(X(i,:),V(i,:));
        c = (norm(X(i,:)))^2-P.R^2;

        dt_starFront = (-b + sqrt(b^2-4*a*c))/(2*a);     % time of impact
        dt_starBack = (-b - sqrt(b^2-4*a*c))/(2*a);
        x_star = X(i,:) + dt_starBack*V(i,:);           % point of impact     

        X_new(i,:) = x_star + (P.dt-dt_starFront)*V(i,:);
    end
end

X = X_new;
V = V_new;
end

