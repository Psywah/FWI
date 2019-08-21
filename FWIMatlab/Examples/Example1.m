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
beta=0;        % Sensor optimisation regularisation parameter 
regp=0;        % regularisation position


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



%% Sensor Optimisation
 
% Starting guess for sensor position  
%p0=model.zr;

% % BB optimisation
%   [pk,hist, error, mw, pvec] = BBOptimOnePixel(p0,1e-12,100,m,model, alpha, m0, beta, regp, fwitol, fwimaxit, mtrue);
 
%  disp(['Optimum sensor location, when x = ' num2str(x(1)) ' is z = ' num2str(pk)])
%  disp(['\psi at this position = ' num2str(hist(end,2))])
 
% %final velocity model
%  vw = reshape(real(1./sqrt(mw)),n);
% %Ground truth model
% vtrue = reshape(real(1./sqrt(mtrue)),n);

% % %Plotting

% %GT
%  figure;
%  imagesc(1e-3*x,1e-3*z,v)
%  axis equal tight
%  title('Ground Truth'); 
%  xlabel('x [km]','fontsize',20);
%  ylabel('z [km]','fontsize',20);


% %Noisy model  
% % % figure;
% % %  imagesc(1e-3*x,1e-3*z,vnoise)
% % %  axis equal tight
% % %  title('Nosiy Model'); 
% % %  xlabel('x [km]','fontsize',20);
% % %  ylabel('z [km]','fontsize',20);

% %Initial 
% %  figure;
% %  imagesc(1e-3*x,1e-3*z,v0)
% %  axis equal tight
% %  title('Initial Velocity');
% %  xlabel('x [km]','fontsize',20);
% %  ylabel('z [km]','fontsize',20);
%   hold on
%   plot(model.xr*1e-3,model.zr*1e-3, 'o-','MarkerFaceColor','red','MarkerEdgeColor','red')
%  plot(model.xs*1e-3,model.zs*1e-3, 'o','MarkerFaceColor','green','MarkerEdgeColor','green')
   
% %Optimised 
%  figure;
%  imagesc(1e-3*x,1e-3*z,vw)
%  axis equal tight
%  title('Optimised Sensor Position'); 
%  xlabel('x [km]','fontsize',20);
%  ylabel('z [km]','fontsize',20);
%  hold on 
% plot(model.xs*1e-3,model.zs*1e-3, 'o','MarkerFaceColor','green','MarkerEdgeColor','green', 'MarkerSize', 8)
%  plot(model.xr*1e-3,p0*1e-3, '*','MarkerFaceColor','red','MarkerEdgeColor','red')
% plot(model.xr*1e-3,pk*1e-3, 'o','MarkerFaceColor','red','MarkerEdgeColor','red', 'MarkerSize', 5) 
% legend('Source', 'Initial Sensor Positions','Optimised Sensor Position', 'Location', 'EastOutside')
%  
 
% %Upper misfit  
% figure
% plot(hist(:,1), hist(:,2)/hist(1,2), 'linewidth', 2)
% xlabel('Iteration Number', 'fontsize', 14)
% ylabel('f(x)=1/2||m_{true}-m_w||^2', 'fontsize', 14)
% title('Function Value', 'fontsize', 16)
% 

% %Upper gradient
% figure
% plot(hist(:,1), hist(:,3)/hist(1,3), 'linewidth', 2)
% xlabel('Iteration Number', 'fontsize', 14)
% ylabel('||g(x)||_2', 'fontsize', 14)
% title('Norm of Gradient', 'fontsize', 16)
