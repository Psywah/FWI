function [mw, hist, error]=Continuation(m, model,p, alpha, m0, tol, maxit,type, ngroup, overlap)

%% Frequency continuation 

% m      = ground truth model (used to calculate synthetic data)
% m0     = Initial guess  
% tol    =  fwi tolerance
% maxit  = fwi max iteraions
% alpha  = fwi regularisation 
% p      = sensor positions
% type   = type of frequency continuation 
% ngroup = number of frequency groups
% overlap = size of overlap between groups


%% Type 1 = Separate groups with overlaps (sequential)

if type==1    
   %disp('Sequential')
    
   % Form frequency groups 
   s=round(size(model.f,2)/ngroup)+overlap;        % size of groups
   fgroup=cell(ngroup, 1);                         % cell to store all frequency groups
       for  j=1:ngroup                             % loop for filling cell array 
           if j*(s-overlap)+overlap>size(model.f,2)
               fgroup{j}=model.f(j*(s-overlap)-(s-overlap)+1:end);    
           else
            fgroup{j}=model.f(j*(s-overlap)-(s-overlap)+1:j*(s-overlap)+overlap);
           end
       end
   
    
    % loop over each group
    for k=1:ngroup
        
        % Get frequencies in group
        f=fgroup{k};   
        model.f=f;
        
       
        % FWI for group
        [mw, hist, error]=FWIOnePixel(m, model,p, alpha , m0, tol, maxit); 
  
        % Initial guess for next group is result of previous
        m0=mw;
       
       
       %disp(['group ' num2str(k)])
       
       
    end
       
       
       
       
end

%% Type 2 = Progressively larger groups - (progressive)
%           Add higher frequencies to the same group at each step

 if type==2
    % disp('Progressive')
    
    % Form frequency groups
    s=round(size(model.f,2)/ngroup);     % size of first group
    fgroup=cell(ngroup, 1);              % cell to store freq groups
    
    for j=1:ngroup
        if j*s>size(model.f,2)
               fgroup{j}=model.f(1:end);
        else
       fgroup{j}=model.f(1:j*s); 
        end
    end
    
    
    % loop over each group
    for k=1:ngroup
        
       % Get frequencies in group 
       f=fgroup{k};   
       model.f=f;
       
       % FWI for group
       [mw, hist, error]=FWIOnePixel(m, model,p, alpha, m0, tol, maxit); 
       
       % Initial guess for next group = result of previous
       m0=mw;
       
       %disp(['group ' num2str(k)])
              
       
    end
    
end   





end 