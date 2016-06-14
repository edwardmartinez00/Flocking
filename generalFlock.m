%%
%% generalFlock
%% Simulates 3 zones model of alignment, attraction & repulsion with N birds
%% With boundary conditions for square or circle

%clear;
close;
if ~exist('s','var')
    s = rng;
end
rng(s);
% addpath('/Users/student/Documents/MATLAB/Flocking/BC');
addpath('C:\Users\thatmathguy\Documents\MCTP\MATLAB\Flocking\BC');
addpath('C:\Users\thatmathguy\Documents\MCTP\MATLAB\Flocking\Save');
tic;

%% ---------  Parameters  --------- %%

% Simulation
P.N = 100;          % number of birds
t = 500;            % time to run
P.dt = .1;          % time step

% Initialization 
X0scale = 7;        % scale of variation in X0 (should be less than L and sqrt(2)*R/2)
V0scale = 1;        % scale of variation in V0
Xshift = 1;         % shift of X ([0,2] with no shift, scale 1)
Vshift = 1;         % shift of V

% Bounds
P.R = 10;         % radius of circle boundary
P.L = 10;         % length of apothem

P.absorb = .8;  % absorbtion scale of reflexive BC (percent)
                % absorbancy = 2 - absorb (in function)

% Strength of flocking attributes
P.IOphi = 1;      % toggles(scales) alignment
P.IOpsi_att = 1; % toggles attraction (.1 recommended)
P.IOpsi_rep = 1; % toggles repulsion
P.IOw = 3;        % scales soft boundary repulsion strength
P.IOpro = 1;      % toggles propulsion

% Distances of function effects
P.a = 2;      % repulsion|a|neutral
P.b = 4;      % neutral|b|attraction
P.c = 8;     % attraction|c|neutral
P.d = 5;      % distance for boundary force

% Graph
shouldPlot = 1;     % logical operator, 1 for plotting
save = 0;           % logical operator, 1 for saving into folder and no pause
pauseTime = .01;    % pause time between frames
window = 13;    % size of square window (1+bound recommended)

% Saving
FileLocation='C:\Users\thatmathguy\Documents\MCTP\MATLAB';%where to create new folder
NewFolder='saveframetest';%name of new folder to be created
Title='plot';%name prefix for saved files

%% ------------ Boundary Choice --------------- %%
%{
The first boundary condition with value 1 will be used.
If all values are 0, there will be no boundary conditions.

Impact causes soft-potential repulsion where birds avoid impacting the boundary in front of them.
Reflect causes reflection with absorbtion as specified above.
Soft causes soft-potential repulsion where birds avoid impacting the boundary in front of them 
and the boundary closest to them.
Proximity causes soft-potential repulsion where birds avoid impacting the boundary closest to them.
%}

% circular boundaries       
BCcircImpact = 0;
BCcircReflect = 0;
BCcircSoft = 0;
BCcircProximity = 0;

% square boundaries
BCsqImpact = 0;
BCsqReflect = 1;
BCsqSoft = 0;
BCsqProximity = 0;
%% V' functions

%addpath('/Users/student/Documents/MATLAB/Flocking');
V_functions;

% phi = @(r) 1./(1+r);    % alignment factor
% %psi = @(r,a,b) (r>a).*(r<b);    % atrraction/repulsion factor
% psi = @(r,a,b,c) (-(r<a) + (r>b.*r<c));
%% Initial conditions

X = X0scale*(2*rand(P.N,2)-Xshift);
V = V0scale*(2*rand(P.N,2)-Vshift);

%% Initialize SaveFrame
if exist(strcat(FileLocation,'\',NewFolder),'file')
    rmdir(strcat(FileLocation,'\',NewFolder),'s');  % delete folder if already exists
end
mkdir(FileLocation,NewFolder)%creates new folder

%% Variance allocation
% varX = zeros(t/P.dt+1,1);
% varV = zeros(t/P.dt+1,1);
% meanV = zeros(t/P.dt+1,2);

% pause
%%-----------------------------------------%%
%%--------------- Time loop ---------------%%
for n=0:P.dt:t
    %% Boundary Conditions
    % New Birds
    % Boundaries and system dynamics
    if BCcircImpact
        str = 'circleImpact';       % for creating function pointer
        [X,V] = RK4birds(X,V,P,str);
        theta=linspace(0,2*pi);
        BoundX=P.R*cos(theta);
        BoundY=P.R*sin(theta);
    elseif BCcircReflect
        str = 'noBounds';
        for i = 1:P.N       % Push bird in if outside boundary (Dirty)
            if norm(X(i,:)) > P.R 
                X(i,:) = P.R*X(i,:)/norm(X(i,:));
            end
        end
        [X_new,V_new] = RK4birds(X,V,P,str);
        [X,V] = circleReflex(X,V,X_new,V_new,P);
        theta=linspace(0,2*pi);
        BoundX=P.R*cos(theta);
        BoundY=P.R*sin(theta);
    elseif BCcircSoft
        str = 'circleSoft';
        [X,V] = RK4birds(X,V,P,str);
        theta=linspace(0,2*pi);
        BoundX=P.R*cos(theta);
        BoundY=P.R*sin(theta);
    elseif BCcircProximity
        str = 'circleProximity';
        [X,V] = RK4birds(X,V,P,str);
        theta=linspace(0,2*pi);
        BoundX=P.R*cos(theta);
        BoundY=P.R*sin(theta);
%----------------------------------------------------        
    elseif BCsqImpact
        str = 'squareImpact';
        [X,V] = RK4birds(X,V,P,str);
        BoundX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
        BoundY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)]; 
    elseif BCsqReflect
        str = 'noBounds';
%         for i = 1:P.N
%             if (X(i,1)>P.L)
%                 X(i,1) = P.L;
%             elseif (X(i,1)<-P.L)
%                 X(i,1) = -P.L;
%             elseif (X(i,2)>P.L)
%                 X(i,2) = P.L;
%             elseif (X(i,2)<-P.L)
%                 X(i,2) = -P.L;
%             end
%         end
        [X_new,V_new] = RK4birds(X,V,P,str);
        [X,V] = squareReflex(X,V,X_new,V_new,P);
        BoundX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
        BoundY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)]; 
    elseif BCsqSoft
        str = 'squareSoft';
        [X,V] = RK4birds(X,V,P,str);
        BoundX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
        BoundY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)]; 
    elseif BCsqProximity
        str = 'squareProximity';
        [X,V] = RK4birds(X,V,P,str);
        BoundX=[linspace(-P.L,P.L),-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100)];
        BoundY=[-P.L*ones(1,100),linspace(-P.L,P.L),P.L*ones(1,100),linspace(-P.L,P.L)]; 
    else
        str = 'noBounds';
        [X,V] = RK4birds(X,V,P,str);
        BoundX = [];
    end
    
    %% Variance
%     varX(round(n/P.dt+1)) = varbirds(X);
%     varV(round(n/P.dt+1)) = varbirds(V);
     %meanV(round(n/P.dt+1),:) = mean(V,1)
    
    %% Plotting
    if shouldPlot
        quiver(X(:,1),X(:,2),V(:,1),V(:,2),0)
        axis([-window window -window window],'square');         %window size
        xlabel('x'); ylabel('y');
        title(['t = ', num2str(n,'%10.2f')]);        %title with time stamp
        
        if ~isempty(BoundX)
            hold on;
            plot(BoundX,BoundY,'k')         %plot the boundary
            hold off;
        end
        
        %set(gcf,'position',get(0,'screensize'))     %large figure    
        
        if save
            SaveFrame(FileLocation,NewFolder,Title,round(n/P.dt+1));
            pauseTime = 0;
        end
    
    % Variance in X plot
%     subplot(3,1,2)
%     plot(0:dt:n,varX(1:round(n/dt+1)))                           % ADD 'IF' STATEMENT
%     xlim([0,t]); ylim([0,200]);
%     xlabel('t'); ylabel('var')
%     title('Variance in Position');
    
    % Variance in V plot
%     subplot(3,1,3)
%     plot(0:dt:n,varV(1:round(n/dt+1)))                           % ADD 'IF' STATEMENT
%     xlim([0,t]); ylim([0,2]);
%     xlabel('t'); ylabel('var')
%     title('Variance in Velocity');

    
        shg
        pause(pauseTime)      %pause to make movie
        %drawnow
    end
end

toc;