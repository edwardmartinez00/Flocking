function X_hat = sqboundaryDistance(X,P)
% used in the flocking model with soft potential boundary conditions
% X_hat is an Nx2 matrix of the points on the boundary closest to birds
% X is an Nx2 matrix of the current positions of N birds
% the domain is the square [-L,L]x[-L,L]

X_hat = zeros(P.N,2);

for i = 1:P.N
    if (X(i,1)>=X(i,2)) && (X(i,1)>-X(i,2))     % East wall
        X_hat(i,:) = [P.L X(i,2)];
    elseif (X(i,1)<=-X(i,2)) && (X(i,1)>X(i,2)) % South wall
        X_hat(i,:) = [X(i,1) -P.L];
    elseif (X(i,1)<=X(i,2)) && (X(i,1)<-X(i,2)) % West wall
        X_hat(i,:) = [-P.L X(i,2)];
    else                                        % North wall
        X_hat(i,:) = [X(i,1) P.L];
    end
end
end

