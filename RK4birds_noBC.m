function [X,V] = RK4birds_noBC(X,V,P)
% applies Runge-Kutta (RK4) method to update position and velocity in flocking model
% no boundary conditions
% outputs X and V are updated positions and velocities, respectively
% input X is an Nx2 matrix of current positions
% input V is an Nx2 matrix of current velocities
% all other inputs are parameters explained in FlockingGeneral.m
% strBound is a string to be converted to a function pointer

% addpath('C:\Users\thatmathguy\Documents\MCTP\MATLAB\Flocking\BC');
BC = str2func('noBounds');  % function pointer

if (rem(P.n,1)==0) && (P.n<P.N)    % sparse control
    if (P.Control) && (P.DogExternal)
        P.Vd = ones(P.N,1)*P.Vd;
        P.nu = [ones(P.n,1)*[P.nu P.nu]; ones(P.N-P.n,1)*[1000000 1000000]];
    elseif (~P.Control) && (P.DogInternal)
        P.Control = 1;
        P.nu = 1000000000;
        P.nu = [P.nu_dog P.nu_dog; ones(P.N-1,1)*[P.nu P.nu]];
    elseif (P.Control) && (P.DogInternal)
        P.nu = [P.nu_dog P.nu_dog; ones(P.n-1,1)*[P.nu P.nu]; ones(P.N-P.n,1)*[1000000 1000000]];
    else
        P.Vd = ones(P.N,1)*P.Vd;
        P.nu = [ones(P.n,1)*[P.nu P.nu]; ones(P.N-P.n,1)*[1000000 1000000]];
    end
else                   % control whole flock
    if (P.Control) && (P.DogExternal)
        P.Vd = ones(P.N,1)*P.Vd;
    elseif (~P.Control) && (P.DogInternal)
        P.Control = 1;
        P.nu = 1000000000;
        P.nu = [P.nu_dog P.nu_dog; ones(P.N-1,1)*[P.nu P.nu]];
    elseif (P.Control) && (P.DogInternal)
        P.nu = [P.nu_dog P.nu_dog; ones(P.N-1,1)*[P.nu P.nu]];
    else
        P.Vd = ones(P.N,1)*P.Vd;
    end
end

if (P.DogInternal)
    % calculate distance between dog and each bird
    Dog_X = ones(P.N,1)*X(1,:);
    Mat_Dog = Dog_X-X;
    Dog_D = sqrt((Mat_Dog(:,1)).^2+(Mat_Dog(:,2)).^2);
    
    % ensure that the dog does not follow or repel itself
    M = max(Dog_D);
    Dog_D(1,:) = max(M,P.d)+1;
    
    if (~P.neighbors)
        % find index of closest bird
        [m,i] = min(Dog_D);
        
        % make Vd an Nx2 matrix
        P.Vd = [X(i,1)-X(1,1) X(i,2)-X(1,2); ones(P.N-1,1)*P.Vd];
    else
        % calculate center of mass of birds within neighborhood
        Dog_D(Dog_D>P.neighborhood) = 0;
        X1 = X(:,1);
        X2 = X(:,2);
        centerX = sum(X1(find(Dog_D)))/length(find(Dog_D));
        centerY = sum(X2(find(Dog_D)))/length(find(Dog_D));
        
        % make Vd an Nx2 matrix
        P.Vd = [centerX-X(1,1) centerY-X(1,2); ones(P.N-1,1)*P.Vd];
    end
end

if (P.RK4)
    % first point is current position and velocity
    X1 = X;
    V1 = V;

    % initial slope are current velocity and acceleration
    kX1 = V1;
    kV1 = BC(X1,V1,P);

    % second point halfway through time step
    X2 = X1+P.dt/2*kX1;
    V2 = V1+P.dt/2*kV1;

    % slopes calculated with second point
    kV2 = BC(X2,V2,P);
    kX2 = V1+P.dt/2*kV2;

    % third point halfway through time step
    X3 = X1+P.dt/2*kX2;
    V3 = V1+P.dt/2*kV2;

    % slopes calculated with third point
    kV3 = BC(X3,V3,P);
    kX3 = V1+P.dt/2*kV3;

    % fourth point at end of time step
    X4 = X1+P.dt*kX3;
    V4 = V1+P.dt*kV3;

    % slopes calculated using fourth point
    kV4 = BC(X4,V4,P);
    kX4 = V1+P.dt*kV4;

    % weighted average of slopes (see RK4 method)
    V_RK = (kX1+2*kX2+2*kX3+kX4)/6;
    dV_RK = (kV1+2*kV2+2*kV3+kV4)/6;

    % updated positions and velocities
    X = X+P.dt*V_RK;
    if (P.Control)
        V = (V+P.dt*dV_RK+(P.dt^2).*P.Vd./P.nu)./(1+(P.dt^2)./P.nu);
    else
        V = V+P.dt*dV_RK;
    end
elseif (P.RK2)
    % first point is current position and velocity
    X1 = X;
    V1 = V;

    % initial slope are current velocity and acceleration
    kX1 = V1;
    kV1 = BC(X1,V1,P);

    % second point at end of time step
    X2 = X1+P.dt*kX1;
    V2 = V1+P.dt*kV1;

    % slopes calculated with second point
    kV2 = BC(X2,V2,P);
    kX2 = V1+P.dt*kV2;
    
    % average of slopes (see improved Euler method)
    V_RK = (kX1+kX2)/2;
    dV_RK = (kV1+kV2)/2;

    % updated positions and velocities
    X = X+P.dt*V_RK;
    if (P.Control)
        V = (V+P.dt*dV_RK+(P.dt^2).*P.Vd./P.nu)./(1+(P.dt^2)./P.nu);
    else
        V = V+P.dt*dV_RK;
    end
else
    % Euler method
    V_RK = V;
    dV_RK = BC(X,V,P);
    
    % updated positions and velocities
    X = X+P.dt*V_RK;
    if (P.Control)
        V = (V+P.dt*dV_RK+(P.dt^2).*P.Vd./P.nu)./(1+(P.dt^2)./P.nu);
    else
        V = V+P.dt*dV_RK;
    end
end

