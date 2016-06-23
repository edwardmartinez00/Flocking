%%
%% V' Functions
%% Functions that affect velocity (for use in flocks)

phi = @(r,P) P.IOphi ./ (1+r);         % alignment function
Psi = @(r,P) (-P.IOpsi_rep*(r<P.a)+P.IOpsi_att*(r>P.b).*(r<P.c)); %attraction/repulsion function
w = @(r,P) P.IOw *(r<P.d).*(1./r-1/P.d);            % soft BC function
pro = @(NormSqV,V,P) P.IOpro*V.*(1-NormSqV); %self-propulsion function
w_dog = @(r,P) P.IOw_dog*(r<P.d_dog).*(1./(1+r)-1/(1+P.d_dog)); % dog repulsion function
