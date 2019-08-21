function [S,U, dRu]=Data(m, model,p)

%% Create synthetic data
 
% Inputs:
%   m                  = model [s^2/km^2]
%   model.h            = [dz,dx] gridspacing in z and x direction [m]
%   model.n            = [nz,nx] number of gridpoints in z and x direction
%   model.f            = frequencies
%   model.xs, model.zs = x and z source locations  
%   model.zr/ model.xr = receiver z/x coordinates location  
%   p                  = receiver coordinates location of interest

% Outputs:
%   S   = data as 3D array (Nr x Ns x Nf)
%   U   = wavefield computed with true model  (n x Ns x Nf)
%   dRu = dR/dp*U (Nr x Ns x Nf)

%% Create finer grid avoiding inverse crime

hh=model.h./2;                           % refine grid by 2
nn=model.n.*2-1;

v = reshape(real(1./sqrt(m)),model.n);   % ground truth velocity
vv=interp2(v, 'nearest');                % define velocity on finer grid
mm=1./(vv(:)).^2;                        % model on finer grid

% Original 'inverse crime' code
 hh=model.h; nn=model.n; mm=m;  

%% Compute arrays on fine grid

% Define source as 2D array 
[Q]=LinearInterpQT(model.zs,model.xs,nn, hh); %transmission - constant x

% Define sampling operator R (extracts u at (z,x) receiver positions)
[R,dRdp]=LinearInterpT(p,model.xr,nn, hh);    %transmission - constant x

% Initialise arrays
S = zeros(size(R,1),size(Q,1),length(model.f));
U=zeros(prod(nn),size(Q,1),length(model.f)); 
dRu=zeros(size(dRdp,1), size(Q,1), length(model.f));

% Loop over frequencies fro computations
for k=1:length(model.f)
    
    Ak = getAR(model.f(k),mm,hh,nn);          % Helmholtz matriz
    Uk=Ak\(Q');                               % Wavefield
    D = R*Uk;                                 % Data
    
    S(:,:,k) = D;                             % Save in arrays
    U(:,:,k) = Uk;
    dRu(:,:,k) = dRdp*Uk;

end



end