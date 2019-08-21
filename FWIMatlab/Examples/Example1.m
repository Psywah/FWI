%% FWI example 

close all

%% Define Model
                
h=[25 25]; n=[101 101];         %spacing & grid size 

z=[0:n(1)-1]*h(1);              %z coordinates 0-2.5km
x=[0:n(2)-1]*h(2);              %x coordinates 0-2.5km

v0=2000*ones(n);                %initial velcoity model
dv=zeros(n);                    %velcoity perturbation (obstacle)
kx=31:71; kz=31:71;             %Pixels to perturb 
dv(kz, kx)=100;                 %Velcoity perturbation
v=v0+dv;                        %Ground truth wavespeed
m=1e6./(v0(:)+dv(:)).^2;        %Perturbed slowness model (ground truth)
mtrue=m;                        %Ground truth model
m0=1e6./(v0(:)).^2;             %Initial slowness model (initial guess)

% Add Noise
% % %eta=rand(n);                     % random no.'s from uniform distribution              
% % eps=1e-2*v(1)*eta;                % Noise % of background  
% % vnoise=v+eps;                     % Add noise to velocity
% % m=1e6./(vnoise(:)).^2;            % Noisy model
                  
f=[0.1:0.4:8];                        % frequencies

%% Acquisition 

%Receiver locations 
zr=[0 1250 2500];
xr=2400*ones(size(zr));

% Source locations
zs=[0 1250 2500];
xs=100*ones(size(zs));

%% Set parameters

%model parameters
model.h=h;
model.n=n;
model.f=f;
model.xs=xs;
model.xr=xr;
model.zs=zs;
model.zr=zr;

%FWI parameters 
fwitol=1e-3;
fwimaxit=100;

%% Regularisation
alpha=0;       % FWI regularisation parameter

%% FWI

% %Original FWI   
 t0=tic;
  [mw, hist, error]=FWIOnePixel(m, model,model.zr, alpha, m0, fwitol, fwimaxit);
 toc(t0) 
 
 disp(['value of psi: ', num2str(ObjFn(mw,m,0,0,0))])
 
% % Frequency Continuation FWI
%  
%  ngroup=10;   % # frequency groups
%  overlap=1;  % # elements shared by adjacent group
%  type=2;     % 1 for sequential, 2 for progressive
%         
%  t1=tic;
%  [mw, hist, error]=Continuation(m, model,model.zr, alpha, m0, fwitol, fwimaxit,type, ngroup, overlap);
%  toc(t1)
%  disp(['value of \psi: ', num2str(ObjFn(mw,m,0,0,0))])
%  

 % % final velocity model
  vw = reshape(real(1./sqrt(mw)),n);
  
 % % Plotting 

 %GT
 figure;
 imagesc(1e-3*x,1e-3*z,v)
 axis equal tight
 title('Ground Truth Velocity'); 
 xlabel('x [km]','fontsize',20);
 ylabel('z [km]','fontsize',20);
 
 %Initial
 figure;
 imagesc(1e-3*x,1e-3*z,v0)
 axis equal tight
 title('Initial Velocity');
 xlabel('x [km]','fontsize',20);
 ylabel('z [km]','fontsize',20);

%Reconstruction 
figure;
imagesc(1e-3*x,1e-3*z,vw)
axis equal tight
title('Result');
xlabel('x (km)','fontsize',12);
ylabel('z (km)','fontsize',12);


% figure
% plot(hist(:,1), hist(:,2)/hist(1,2), 'linewidth', 2)
% xlabel('Iteration Number', 'fontsize', 14)
% ylabel('f(x)', 'fontsize', 14)
% title('Function Value', 'fontsize', 16)
% 
% figure
% plot(hist(:,1), hist(:,3)/hist(1,3), 'linewidth', 2)
% xlabel('Iteration Number', 'fontsize', 14)
% ylabel('||g(x)||_2', 'fontsize', 14)
% title('Norm of Gradient', 'fontsize', 16)
% 
%  figure
%  plot(hist(:,1), error, 'LineWidth', 2)
%  xlabel('Iteration Number', 'fontsize', 14)
% ylabel('||m_{true}-m_{w}||', 'fontsize', 14)
% title('Error', 'fontsize', 16)




