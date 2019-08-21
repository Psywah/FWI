function [f,g,H] = Misfit(m,D,alpha,model,p)

%% Evaluate least-squares misfit & gradient & Hessian 
%
%   0.5||RA^{-1}(m)Q - D||_{F}^2 + 0.5\alpha||Lm||_2^2,
%
% where R, Q encode the receiver and source locations and L is the
% first-order FD matrix
%
% input:
%   m             = model - squared-slownes [s^2/km^2]
%   D             =  synthetic data matrix
%   alpha         = regularization parameter
%   model.h       = gridspacing in each direction d = [d1, d2];
%   model.n       = number of gridpoints in each direction n = [n1, n2]
%   model.f       = frequency [Hz].
%   model.zr/xr   = z locations of receivers [m] 
%   p             = receiver location of interest
%   model.{zs,xs} = {z,x} locations of sources [m] 
%
% output:
%   f - value of misfit
%   g - gradient (vector of size size(m))
%   H - GN Hessian (function handle)


%% Get matrices

%first derivative matrix
L = getLR(model.h,model.n);

%Receiver Sampling
%[R,~]=LinearInterp(p,model.zr,model.n, model.h);  %Reflection
[R,~]=LinearInterpT(p,model.xr,model.n, model.h);  %Transmssion

%Source Sampling
%[Q]=LinearInterpQ(model.xs,model.zs,model.n, model.h); %Reflection
[Q]=LinearInterpQT(model.zs,model.xs,model.n, model.h); % Transmission

m = m(:);  %write model as vector


% initialize misfit and gradient
f = alpha*0.5*norm(L*m)^2; % regularisation term
g = alpha*(L'*L)*m;        % differentaite reg term

% lopp over frequencies
for k = 1:length(model.f)
    
	% get forward Helmholtz operator
    Ak = getAR(model.f(k),m,model.h,model.n);
    
    % gradient of forward operator * wavefield
    % G= (dA/dm u)
    Gk = @(u)getGR(model.f(k),m,u,model.h,model.n);
    
	% solve forward wave-equation 
	Uk = Ak\(Q');
    
	% solve adjoint wave-equation [ lam= A*^{-1} R^T(d-Ru) ]
    Vk = Ak'\(R'*(D(:,:,k) - R*Uk));
    
	% compute misfit f (summed over frequencies)
    f = f+0.5*norm(R*Uk - D(:,:,k),'fro')^2;   

    % Compute gradient g =  Re[(dA/dm u)* lambda] = Re[G* lambda]
    
    % Sum over Sources
    for j = 1:size(Uk,2)
        g = g + real(Gk(Uk(:,j))'*Vk(:,j));   
    end
    
    % Gradient is real
    g=real(g);

end


%% get H (GN Hessian)
H = @(x)Hmv(x,m,alpha,model);

end

function y = Hmv(x,m,alpha,model)

% get matrices
L = getLR(model.h,model.n);

[R,~]=LinearInterp(p,model.zr,model.n, model.h);
[Q,~]=LinearInterp(model.xs,model.zs,model.n, model.h);


% compute mat-vec
y = alpha*(L'*L)*x;

% loop over frequencies
for k = 1:length(model.f)
Ak = getAR(model.f(k),m,model.h,model.n);    
Gk = @(u)getGR(model.f(k),m,u,model.h,model.n);
Uk = Ak\(Q');

% loop over sources 
for j = 1:size(Uk,2);
    y = y + real(Gk(Uk(:,j))'*(Ak'\((R'*R)*(Ak\(Gk(Uk(:,j))*x)))));
end

end
end

