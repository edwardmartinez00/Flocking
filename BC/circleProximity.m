%% dV = circleProximity(X,V,P)
%
% Finds change in velocity at a soft circle boundary based solely on the
% closest point to the birds
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position depending on distance to boundary
% P.N, number of birds
% P.R, radius of circle
% P.dt, time step size
function dV = circleProximity(X,V,P)
    V_functions;

    dV = noBounds(X,V,P);
    
    % allocating space
    X_hat = zeros(P.N,2);
    
    for i = 1:P.N     % Bird loop
        X_hat(i,:) = P.R*X(i,:)/norm(X(i,:));        % closest wall
    end   
    
    % Calculating distance to closest wall
    D_hat = X_hat-X;
    r_hat = sqrt(sum(D_hat.^2,2));
    r_hat = r_hat*[1 1];
    
    dV = dV - w(r_hat,P).*(D_hat./r_hat);
    
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