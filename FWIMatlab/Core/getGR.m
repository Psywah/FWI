function G = getGR(f,m,u,h,n)

%% Define matrix G(m,u) = d(A(m))/dm*u

% Input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   u - wavefield
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   G - sparse matrix

% Angular frequency
omega = 1e-3*2*pi*f;

% No. of discretisation points 
N     = prod(n);

a = ones(n); a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

% Define sparse matrix
G = omega^2*diags(a.*u) + (2i*omega/h(1))*diags(0.5*(1-a).*u./sqrt(m));

end
