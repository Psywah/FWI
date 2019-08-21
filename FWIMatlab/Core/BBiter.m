function [xk,hist] = BBiter(fh,x0,tol,maxit)

%% Optimisation for FWI

% Barzalai Borwein iteration to solve 
%       min_m ||d-Ru(m)||  
 
%  Input:
%   fh    = function handle that returns value and gradient: [f,g] = fh(x)
%   x0    = initial iterate / guess
%   tol   = stop when ||g||_2 <= tol
%   maxit = stop when iter > maxit
%   m     = ground truth model
 
%  Output:
%   xk    = final iterate of model
%   hist  = array with rows [iter, f, g]


% Initial guess 
k       = 0;
xk      = x0;
[fk,gk] = fh(xk);
cnst=10^5;
tk=norm(gk)*cnst; 


% Save history
hist = [k,fk,norm(gk)];
hist0=hist;
 

% display initial f(x0) & g(x0)
% %%%%%%%%%%%%%%
% fprintf('----------------------------------- \n')
% 
% fprintf(1,' k , fk          , ||gk||_2\n');
% fprintf(1,'%3d, %1.5e , %1.5e \n',hist);
%%%%%%%%%%%%%%

 

% update model until norm(gradient) >= tolerance or max iterations exceeded   
while (norm(gk) > tol)&&(k < maxit)  
         
     sk=-gk/tk;     % step 
     xk = xk + sk;  % update model
     
     
     % gradient evaluation
     [fk,gk] = fh(xk);
   
     % update steplength
     tk=tk+(sk'*gk)/norm(sk)^2;             
   
   
     % update counter
     k = k + 1;
   
     % update history
     hist(k,:) = [k,fk,norm(gk)];
     
     %%%%%%%%
     % display history 
     % fprintf(1,'%3d, %1.5e , %1.5e \n',k,fk,norm(gk));
     %%%%%%%%
   
end

%fprintf('----------------------------------- \n')

% output history
hist=[hist0; hist];


end



 % update BB step and model pixel
%    sk = -gk(ind)/tk;
%    xk(ind) = xk(ind) + sk;
     
%%one value of g
     %sk=-gk(ind(1))/tk; %one value of gradient
     %xk(ind) = xk(ind) + sk;  