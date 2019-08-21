function [R,dRdp]=LinearInterpT(p,xr,N, h)

%% Linear Interpolant in z direction

% Crosswell/ transmission type experiments (constant in x)

% Function file that intakes sensor position for multiple sensors 
% and outputs the interpolant matrix (interpolating in x and z direction)
% Assumes that x coord of sensors is constant (all along a line in x)

% Inputs: 
%   p  = z sensor position, 
%   xr = x sensor position, 
%   no. nodes N, spzcing h

% Outputs: 
% R    = interpolant matrix, 
% dRdp = derivative of interpolant derivative 

% Set up for Transmission experiments - varied z and constant x
% For reflection experiments swap order of inputs

%% Interpolate in z direction

    Nr=length(p);           % nummber of receivers 
    
    z  = [0:N(1)-1]*h(1);   % spatial vectors
    x  = [0:N(2)-1]*h(2);
    
    Iz=zeros(Nr,N(1));      % Initialise z interpolant (Nr x Nx)
    dIzdp=zeros(Nr,N(1));   % Initialise derivative of z interpolant
    
    % Loop over sensors 
for i=1:Nr

        % find the index of the closest node iz
        [~,iz]=min(abs(z-p(i))); 
        
        % find value of z at that node
        z1=z(iz);


        % find the node on the other side of p
        if z1>p(i)            % if the known node is on the right of p
            z2=z1;            % rename right node and index
            iz2=iz;
            iz1=iz-1;         % define left node and index
            z1=z(iz1);
        elseif z1<p(i)        % else if known node is on left
            iz2=iz+1;         % define right node and index
            z2=z(iz2);
            iz1=iz;           % define left index
        elseif z1==p(i)       % if receiver is on node (discontinuity problem) 
            iz1=iz;            
            Iz(i,iz1)=1;      % Interp matrix is value one at node    
            if z1==z(end)     % if p is at last node on right
                z2=z1;        % define dRdp as left derivative  
                iz2=iz1;
                iz1=iz2-1;
                z1=z(iz1);
            else              % else define as right deriavtive
                iz2=iz+1;
                z2=z(iz2);
            end
          dIzdp(i,iz1)=-1/(z2-z1);
          dIzdp(i,iz2)=1/(z2-z1);
            continue
        end

        % Linear interpolant matrix
        Iz(i,iz1)=(z2-p(i))/(z2-z1);
        Iz(i,iz2)=(p(i)-z1)/(z2-z1);

        % Derivative of linear interpolant matrix
        dIzdp(i,iz1)=-1/(z2-z1);
        dIzdp(i,iz2)=1/(z2-z1);
        

     
end
    

%% Interpolate in x direction

% Assume x is constant, i.e. x(1)=x(2)=...=x(end) 

    Ix=zeros(1,N(2));      % initialise x-interpolant
   
        % find the index of the closest x node ix
        [~,ix]=min(abs(x-xr(1))); 
        
        % find value of x at that node
        x1=x(ix);
        

        % find the node on the other side of p
        if x1>xr(1)            % if the known node is below zr
            x2=x1;             % rename node below and index
            ix2=ix;
            ix1=ix-1;          % define node above and index
            x1=x(ix-1); 
        elseif x1<xr(1)        % else if known node is above
            ix2=ix+1;          % define right node and index
            x2=x(ix2);
            ix1=ix;            % define above index
        elseif x1==xr(1)       % if receiver is on node (discnty problem) 
            ix1=ix;            
            Ix(1,ix1)=1;       % Interp matrix is value one at node    
            if x1==x(end)     % if zr is at last node at the bottom
                x2=x1;        % define dRdp as derivative above  
                ix2=ix1;
                ix1=ix2-1;
                x1=x(ix1);
            else              % else define as deriavtive below
                ix2=ix+1;
                x2=x(ix2);
            end
        end

        % Linear interpolant matrix
        Ix(1,ix1)=(x2-xr(1))/(x2-x1);
        Ix(1,ix2)=(xr(1)-x1)/(x2-x1);
        

        % Derivative of linear interpolant matrix with respect to x
        % coordinate is zero
        
        
    R=kron(Ix, Iz);                % define R as x & z interpolant
    R=full(R);                     % Size Nr x N^2
     
    dRdp=kron(Ix, dIzdp);         
    dRdp=full(dRdp);
    
end