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
    Mat_xdir = Mat_dXx./(Mat_r+diag(ones(P.N,1)));
    Mat_ydir = Mat_dXy./(Mat_r+diag(ones(P.N,1)));
    
    % Self-propulsion
    NormSqV = V(:,1).^2+V(:,2).^2;  % Norm squared of V
%     Mat_NormSqV = [NormSq_dV,NormSq_dV];
    proVx = pro(NormSqV,V(:,1),P);
    proVy = pro(NormSqV,V(:,2),P);
    
    % Change in velocity calculation
    dVx = 1/P.N*sum(Mat_phi.*Mat_dVx + Mat_psi.*Mat_xdir)'+proVx;
    dVy = 1/P.N*sum(Mat_phi.*Mat_dVy + Mat_psi.*Mat_ydir)'+proVy;
    dV = [dVx, dVy];  % change in V
end    