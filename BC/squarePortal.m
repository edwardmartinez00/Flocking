function [X,V] = squarePortal(X,V,X_new,V_new,P)
% teleports birds
% Detailed explanation goes here

for i = 1:P.N
    if X_new(i,1)>P.L
        X_new(i,1) = X_new(i,1)-2*P.L;
    end
    if X_new(i,1)<-P.L
        X_new(i,1) = X_new(i,1)+2*P.L;
    end
    if X_new(i,2)>P.L
        X_new(i,2) = X_new(i,2)-2*P.L;
    end
    if X_new(i,2)<-P.L
        X_new(i,2) = X_new(i,2)+2*P.L;
    end
end
    
X = X_new;
V = V_new;
end

