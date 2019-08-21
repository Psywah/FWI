function [mw, hist]=FWIfull(m, model,p, alpha, m0, tol, maxit)

% File which carries out all FWI steps 

% Inputs: 
%   m     = ground truth model
%   m0    = initial guess of model 
%   p     = sensor position (1 coordinate)
%   alpha = regularisation parameter
%   tol   = tolerance (stopping criteria)
%   maxit = maximum iteration before stopping optimisation
%   
% Outputs: 
%   mw   = FWI output model / reconstruction   
%   hist = history of misfit and gradient per iteration
%

%% Create synthetic data 
[D, ~, ~] = Data(m,model,p);

%% Find misfit d-Ru (synthetic data - modelled data)        
fh = @(m)Misfit(m,D,alpha,model,p);

%% Optimisation
%disp('FWI Iterations')
[mw,hist] = BBiter(fh,m0,tol,maxit);


end