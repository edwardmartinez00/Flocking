%% dV = squareImpact(X,V,P)
%
% Finds change in velocity at a soft square boundary solely based on impact 
% as a reference
%
% X, position vector
% V, velocity vector
% P.N, number of birds
% P.L, apothem of square
% P.d, distance that W occurs in
% P.dt, time step size
function dV = squareImpact(X,V,P)
    V_functions;

    % calculate acceleration without boundary conditions
    dV = noBounds(X,V,P);
    
    % point on boundary in front of birds
    X_star = sqboundaryImpact(X,V,P);        % would be impact point
    
    % Calculating distance to impact
    D_star = X_star-X;
    r_star = sqrt(sum(D_star.^2,2));
    r_star = r_star*[1 1];    
    
    dV = dV - w(r_star,P).*(D_star./r_star);
    
    % calculate distances between dog and each bird for repulsion purposes
    if (P.DogExternal)
        % calculate distance between dog and each bird
        Dog_X = ones(P.N,1)*P.X_dog;
        Mat_Dog = Dog_X-X;
        Dog_D = sqrt((Mat_Dog(:,1)).^2+(Mat_Dog(:,2)).^2);
        
        % ensure that the dog does not follow or repel itself
        M = max(Dog_D);
        Dog_D(1,:) = max(M,P.d)+1;
    elseif (P.DogInternal)
        % calculate distance between dog and each bird
        Dog_X = ones(P.N,1)*X(1,:);
        Mat_Dog = Dog_X-X;
        Dog_D = sqrt((Mat_Dog(:,1)).^2+(Mat_Dog(:,2)).^2);
        
        % ensure that the dog does not follow or repel itself
        M = max(Dog_D);
        Dog_D(1,:) = max(M,P.d)+1;
    end
    
    % if there is a dog in the model, repel the birds
    if (P.Dog)
        Dog_D = Dog_D*ones(1,2);
        dV = dV-P.sD.*w_dog(Dog_D,P).*(Mat_Dog./Dog_D);
    end
end