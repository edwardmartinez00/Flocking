%% dV = noBounds(X,V,P)
% Finding change in V without bounds
% Set up matrices of differences (xj-xi) and (vj-vi) in components
    
function dV = noBounds(X,V,P)
    V_functions;
    
    onecol = ones(size(X,1),1);
    Mat_dXx = X(:,1)*onecol'-onecol*(X(:,1))'; % matrices of diffs in X
    Mat_dXy = X(:,2)*onecol'-onecol*(X(:,2))';
    Mat_dVx = V(:,1)*onecol'-onecol*(V(:,1))'; % matrices of diffs in V
    Mat_dVy = V(:,2)*onecol'-onecol*(V(:,2))';
    
    % Setting up matrix of phi and psi
    Mat_r = sqrt(Mat_dXx.^2 + Mat_dXy.^2);    % norm of (xj-xi)
    
    Mat_phi = phi(Mat_r,P);
    Mat_psi = Psi(Mat_r,P);
    
    % Direction of attraction/repulsion
    Mat_xdir = Mat_dXx./(Mat_r+diag(ones(P.N,1))+realmin);
    Mat_ydir = Mat_dXy./(Mat_r+diag(ones(P.N,1))+realmin);
    
    % Self-propulsion
    NormSqV = V(:,1).^2+V(:,2).^2;  % Norm squared of V
%     Mat_NormSqV = [NormSq_dV,NormSq_dV];
    proVx = pro(NormSqV,V(:,1),P);
    proVy = pro(NormSqV,V(:,2),P);
    
    % Change in velocity calculation
    dVx = 1/P.N*sum(Mat_phi.*Mat_dVx + Mat_psi.*Mat_xdir)'+proVx;
    dVy = 1/P.N*sum(Mat_phi.*Mat_dVy + Mat_psi.*Mat_ydir)'+proVy;
    dV = [dVx, dVy];  % change in V
    
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