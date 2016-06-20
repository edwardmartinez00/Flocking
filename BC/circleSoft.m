%% dV = circleSoft(X,V,P)
%
% Finds change in velocity at a soft circle boundary
%
% X, position vector (2 columns)
% V, velocity vector (2 columns)
% X_new, potential position depending on distance to boundary
% P.N, number of birds
% P.R, radius of circle
% P.dt, time step size
function dV = circleSoft(X,V,P)
    V_functions;

    dV = noBounds(X,V,P);
    
    % allocating space
    X_star = zeros(P.N,2);
    X_hat = zeros(P.N,2);
    
    for i = 1:P.N     % Bird loop
        a = norm(V(i,:))^2;
        b = 2*dot(X(i,:),V(i,:));
        c = (norm(X(i,:)))^2-P.R^2;

        dt_star = (-b+sqrt(b^2-4*a*c))/(2*a);        % impact time
        X_star(i,:) = X(i,:) + dt_star*V(i,:);       % impact point
        X_hat(i,:) = P.R*X(i,:)/norm(X(i,:));        % closest wall

%         hatdiff = X_hat-X(i,:);
%         r = norm(X(i,:)-X_star);            % distance to impact point

%         X(i,:) = x_star + (P.dt-dt_star)*V(i,:);
%         V(i,:) = V(i,:) - P.dt*w(r,P)*(hatdiff)/norm(hatdiff);
%         dV(i,:) = flock_dV(i,:) - w(r,P)*(hatdiff)/norm(hatdiff);
    end  
    
    % Calculating distance to impact
    D_star = X_star-X;
    r_star = sqrt(sum(D_star.^2,2));
    r_star = r_star*[1 1];    
    
    % Calculating distance to closest wall
    D_hat = X_hat-X;
    r_hat = sqrt(sum(D_hat.^2,2));
    r_hat = r_hat*[1 1];
    
    dV = dV - w(r_star,P).*(D_hat./r_hat);
    
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