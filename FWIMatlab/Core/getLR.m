function L = getLR(h,n)

%% Define 2D first-order FD matrix 

% use:
%   L = getL(h,n)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%
% output
%   L     - sparse first derivative matrix [Dz; Dx] 

    D1 = spdiags(ones(n(1),1)*[-1 1]/h(1),[0:1],n(1)-1,n(1)); %Dz
    D2 = spdiags(ones(n(2),1)*[-1 1]/h(2),[0:1],n(2)-1,n(2)); %Dx
    
    L  = [kron(speye(n(2)),D1); kron(D2,speye(n(1)))];        %[Dz; Dx];
    
end