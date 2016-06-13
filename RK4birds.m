function [X,V] = RK4birds(X,V,P,strBound)
% applies Runge-Kutta (RK4) method to update position and velocity in flocking model
% no boundary conditions
% outputs X and V are updated positions and velocities, respectively
% input X is an Nx2 matrix of current positions
% input V is an Nx2 matrix of current velocities
% all other inputs are parameters explained in FlockingGeneral.m
% strBound is a string to be converted to a function pointer

addpath('C:\Users\thatmathguy\Documents\MCTP\MATLAB\Flocking\BC');
BC = str2func(strBound);  % function pointer

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
V = V+P.dt*dV_RK;
end