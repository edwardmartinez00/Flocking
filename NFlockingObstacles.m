%NflockingGeneral
%flocking model with N birds
%with updated phi alignment and psi attraction/repulsion
%with boundary conditions with 4 obstacles
tic
%Parameters of simulation ---------------------------------------------
N=200; %number of birds
R=10; %radius of boundary circle
L=10; %half side length of boundary square

%------Obstacles----------------------------
%ObstaclesIC=[9*ones(19,1);8*ones(19,1);7*ones(19,1);6*ones(19,1);5*ones(19,1);4*ones(19,1);3*ones(19,1);2*ones(19,1);1*ones(19,1),(-9:9)';(-9:9)';(-9:9)';(-9:9)';(-9:9)';(-9:9)';(-9:9)';(-9:9)';(-9:9)',.25*ones(171,1)];
length=1.5;
%ObstaclesIC=[-5 -5 length; -5 0 length; -5 5 length; 0 -5 length; 0 0 length; 0 5 length; 5 -5 length;
    %5 0 length; 5 5 length];
%ObstaclesIC=[-6 -6 length; -6 -2 length; -6 2 length; -6 6 length; -2 -6 length; -2 -2 length; -2 2 length; -2 6 length; 2 -6 length; 2 -2 length; 2 2 length; 2 6 length; 6 -6 length; 6 -2 length; 6 2 length; 6 6 length];   
%ObstaclesIC=[-6 -6 2*rand; -6 -2 2*rand; -6 2 2*rand; -6 6 2*rand; -2 -6 2*rand; -2 -2 2*rand; -2 2 2*rand; -2 6 2*rand; 2 -6 2*rand; 2 -2 2*rand; 2 2 2*rand; 2 6 2*rand; 6 -6 2*rand; 6 -2 2*rand; 6 2 2*rand; 6 6 2*rand]; 
%ObstaclesIC=[7*(2*rand(16,2)-1),2*rand(16,1)];
ObstaclesIC=[linspace(-10,10)' linspace(-10,10)' ones(100,1)*.1];
L1=1;%half side length of first obstacle
C1=[-5,3];%center of first obstacle
L2=2;%half side length of second obstacle
C2=[-2,-3];%center of second obstacle
L3=3;%half side length of third obstacle
C3=[3,-5];%center of third obstacle
L4=.5;%half side length of fourth obstacle
C4=[3,0];%center of fourth obstacle
    obstacles=size(ObstaclesIC);
    obstacles=obstacles(1);

    
%---------Bird Initial Conditions--------------------------%    
XVariation=10; %amplitude of variation in X(0)
if XVariation>L
    XVariation=L;
end
VVariation=1; %amplitude of variation in V(0)
Window=11; %X and Y max for window
max_time=500; %how long to run
speed=1; %speed of simulation, x0.01 sec
dt=.1; %step size

Q=1; %with alignment,attraction=1, otherwise =0
P=1; %self-propulsion factor
aa=1;%transition repulsion->neutral
bb=5;%transition neutral->attraction
d=3; %distance to start feeling boundary
absorbancy=2; % 1<=absorbancy<2 Absorbancy of BC
sigma=0; %magnitude of noise


%-------V' component functions------------------------------
phi=@(r) Q*10./(1+r); %alignment decay function
psi=@(r) Q*.5*(-1./(r+1)*(r<aa)+(r>aa).*(r<bb)+(r>bb)); %attraction/repulsion function
w=@(r) 1./r*(r<d); %soft BC function
p=@(NormSqV,v) v.*(1-NormSqV); %self-propulsion function
%----------------------------------------------------------------------

%-------------------Initial Conditions-------------------------------------


%X=-1*XVariation+2*XVariation*rand(N,2);%uniform square distribution of position
%X=PosDist(N,XVariation,C1,C2,C3,C4,L,L1,L2,L3,L4);

X=NPosDist(N,XVariation,ObstaclesIC);
%X=CheckDist(N,X,XVariation,ObstaclesIC);
rng(3);
V=-1*VVariation+2*VVariation*rand(N,2);%uniform square distribution of velocity

%-------------------------------RECORDING FRAMES---------------------------
record=0;%1 if record, 0 if not
FileLocation='/Volumes/Samsung USB/MCTP'; %where to create new folder
%FileLocation='D:\MCTP';
NewFolder='Obstacles Video 1'; %name of new folder to be created
Title='plot'; %prefix of filenames
if record==1 
    mkdir(FileLocation,NewFolder)%creates new folder
end
%-------------------------------------------------------------------------

close
close
shg
addpath ('/Volumes/Samsung USB/MCTP/Boundary Conditions')%for mac
addpath ('D:\MCTP/Boundary Conditions')%for Lenovo
addpath ('\\Client\D$\MCTP\Boundary Conditions')%for Asus
addpath ('C:\Users\Kevin\Desktop\Flocking Model\Boundary Conditions')%for Lenovo desktop
spin=0;%initial value for calculating spin
noise=zeros(N,2);%initial value for noise
DT_STAR=[];

BoundaryX=[linspace(-L,L,1000),-L*ones(1,1000),linspace(-L,L,1000),L*ones(1,1000)];
BoundaryY=[-L*ones(1,1000),linspace(-L,L,1000),L*ones(1,1000),linspace(-L,L,1000)];
% obstacle1=square(L1,C1);
% obstacle2=square(L2,C2);
% obstacle3=square(L3,C3);
% obstacle4=square(L4,C4);
 
 %OBSTACLES=[obstacle1 obstacle2 obstacle3 obstacle4];
 OBSTACLES=[];
 for n=1:obstacles
 OBSTACLES=[OBSTACLES, square(ObstaclesIC(n,3),ObstaclesIC(n,[1 2]))];
 end
 

for kstep=1:max_time/dt %overall time loop
   %X=CheckDist(N,X,XVariation,ObstaclesIC);

    %--------noise--------------------------------------------------%
    noise=noise+sqrt(dt)*randn(N,2);
    %=sqrt(dt)*cumsum(rand(1,N));
    
    
    %----------------------------Calculations to update V------------------------------------------------------------%
    mat_r=sqrt((X(:,1)*ones(1,N)-ones(N,1)*X(:,1)').^2+(X(:,2)*ones(1,N)-(ones(N,1)*X(:,2)')).^2);%magnitude of difference in position between i and all other birds
    mat_pos_difx=(X(:,1)*ones(1,N)-ones(N,1)*X(:,1)');%dif between x components of pos
    mat_pos_dify=(X(:,2)*ones(1,N)-ones(N,1)*X(:,2)');%dif between y components of pos
    mat_dirx=mat_pos_difx./(mat_r+diag(ones(N,1)));%direction to align x component
    mat_diry=mat_pos_dify./(mat_r+diag(ones(N,1)));%direction to align y component
  
       
   %---------for alignment-----------------------%
    mat_phi=phi(mat_r);%apply alignment function to distances
    mat_dvx=V(:,1)*ones(1,N)-ones(N,1)*V(:,1)';%dif between x components of velocity        
    mat_dvy=V(:,2)*ones(1,N)-ones(N,1)*V(:,2)';%dif between y components of velocity  
    
    dvx_align=1/N*sum(mat_phi.*mat_dvx);%calculate x component of differential alignment
    dvy_align=1/N*sum(mat_phi.*mat_dvy);%calculate y component of differential alignment
    
    %----------for attraction/repulsion--------------%
    mat_psi=psi(mat_r);%apply attraction/repulsion function to distances
    dvx_attract=1/N*(sum(mat_psi.*mat_dirx));%calculate x component of differential attraction/repulsion
    dvy_attract=1/N*(sum(mat_psi.*mat_diry));%calculate y component of differential attraction/repulsion
    
    %-------for self-propulsion---------------------%
    mag_v_sq=(V(:,1).^2+V(:,2).^2);%magnitude of velocities squared
    mag_v_sq=[mag_v_sq, mag_v_sq];%creates Nx2 matrix ov mag_v_sq
    mat_p=p(mag_v_sq,V);%appliy self propulsion function
    
    
    
    
   %------------------total differential components--------------------------------------% 
    dv=[dvx_align',dvy_align']+[dvx_attract',dvy_attract']+P*mat_p;%+sigma*noise; %total differential components of V
    
   V=V+dt*dv+sigma*noise;%update velocity using Euler method and add noise
  
    %Spin%execute spin calculatio script
  %------------------BC Warning--------------------------%
    Xnew=X+dt*(V);%+sigma*noise); %what next X would be (for boundary warning)
    
    %---------SELECT A BOUNDARY CONDITION------------------------------%
  
    %BCsquareImpactOneObstacle
  %BCsquareImpactObstacles
%BCsquareAbsorbantObstacles
BCsquareAbsorbantNObstacles
  
  




   %-------------Plots----------------------------------%
    plot(BoundaryX,BoundaryY,'.k',OBSTACLES(1,:),OBSTACLES(2,:),'.k',X(:,1),X(:,2),'.')%plot the boundary
    hold on
    %quiver(X(:,1),X(:,2),.5*V(:,1),.5*V(:,2),0,'r')%plot the birds
    title(['Time t=' num2str(kstep*dt,'%10.2f')])%title the plot with time
    set(gca,'FontSize',25)
    axis([-1*Window Window -1*Window Window],'square')
    axis off

    %set(gcf,'position',get(0,'screensize'))%opens plot in full screen
    
    if record==1
        SaveFrame
    end
    %pause(.01*speed)
    drawnow%update plot immediately
    hold off
end
toc

            