%% [X_new,V_new] = squareReflexFail(X,V,X_new,P)
%
% Checks points at a reflexive square boundary
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position vector (depending on relation to bounds)
% P.L, apothem of square
% P.dt, time step size
%
function [X_new,V_new,BoundaryX,BoundaryY] = squareReflexFail(X,V,X_new,P)
    V_new = V;
    absorbancy = 2 - P.absorb;
    
    % Logical indeces of exceeding bounds
    right = X_new(:,1) > P.L;
    left = X_new(:,1) < -P.L;
    up = X_new(:,2) > P.L;
    down = X_new(:,2) < -P.L;
    
    dt_star = zeros(size(X,1),1);
    X_star = zeros(size(X,1),2);
    
    zerocol = zeros(size(X,1),1);
    
    if any(right)
        dt_star(right) = (P.L-X(right,1))./V(right,1);
        mat_dt_star = [dt_star(right),dt_star(right)];
        X_star(right,:) = X(right,:) + (mat_dt_star).*V(right,:);
        V_new(right,:) = V(right,:) - absorbancy*[V(right,1),zerocol(right)];
        X_new(right,:) = X_star(right,:) + (P.dt - mat_dt_star).*V_new(right,:);
    
    elseif any(left)
        dt_star(left) = (-P.L-X(left,1))./V(left,1);
        mat_dt_star = [dt_star(left),dt_star(left)];
        X_star(left,:) = X(left,:) + (mat_dt_star).*V(left,:);
        V_new(left,:) = V(left,:) - absorbancy*[V(left,1),zerocol(left)];
        X_new(left,:) = X_star(left,:) + (P.dt - mat_dt_star).*V_new(left,:);
  
    elseif any(up)
        dt_star(up) = (P.L-X(up,2))./V(up,2);
        mat_dt_star = [dt_star(up), dt_star(up)];
        X_star(up,:) = X(up,:) + (mat_dt_star).*V(up,:);
        V_new(up,:) = V(up,:) - absorbancy*[zerocol(up),V(up,2)];
        X_new(up,:) = X_star(up,:) + (P.dt - mat_dt_star).*V_new(up,:);

    elseif any(down)
        dt_star(down) = (-P.L-X(down,2))./V(down,2);
        mat_dt_star = [dt_star(down),dt_star(down)];
        X_star(down,:) = X(down,:) + (mat_dt_star).*V(down,:);
        V_new(down,:) = V(down,:) - absorbancy*[zerocol(down),V(down,2)];
        X_new(down,:) = X_star(down,:) + (P.dt - mat_dt_star).*V_new(down,:);
    end  
    BoundaryX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
    BoundaryY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)];
 end
