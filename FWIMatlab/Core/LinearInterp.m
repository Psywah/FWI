function [R,dRdp]=LinearInterp(p,zr,N, h)

%% Linear Interpolant in x direction

% Function file that intakes sensor position for multiple sensors 
% and outputs the interpolant matrix (interpolating in x and z direction)
% Assumes that z corrdinate of the sensor is constant

% Inputs: 
%   p  = x sensor position, 
%   zr = z sensor position, 
%   no. nodes N, spzcing h

% Outputs: R interpolant matrix, dRdp derivative of interpolant derivative 

% Set up for reflection experiments - varied x and constant z
% For transmission experiments swap order of inputs

%% Interpolate in x direction

    Nr=length(p);           % nummber of receivers 
    
    z  = [0:N(1)-1]*h(1);   % spatial vectors
    x  = [0:N(2)-1]*h(2);
    
    Ix=zeros(Nr,N(2));      %initialise x interpolant (Nr x Nx)
    dIxdp=zeros(Nr,N(2));
    
    % loop over sensors 
for i=1:Nr

        [~,ix]=min(abs(x-p(i))); 
        % find the index of the closest node ix

        x1=x(ix);
        % find value of x at that node

        % find the node on the other side of p
        if x1>p(i)            % if the known node is on the right of p
            x2=x1;            % rename right node and index
            ix2=ix;
            ix1=ix-1;         % define left node and index
            x1=x(ix1);
        elseif x1<p(i)        % else if known node is on left
            ix2=ix+1;         % define right node and index
            x2=x(ix2);
            ix1=ix;           % define left index
        elseif x1==p(i)       % if receiver is on node (discnty problem) 
            ix1=ix;            
            Ix(i,ix1)=1;      % Interp matrix is value one at node    
            if x1==x(end)     % if p is at last node on right
                x2=x1;        % define dRdp as left derivative  
                ix2=ix1;
                ix1=ix2-1;
                x1=x(ix1);
            else              % else define as right deriavtive
                ix2=ix+1;
                x2=x(ix2);
            end
          dIxdp(i,ix1)=-1/(x2-x1);
          dIxdp(i,ix2)=1/(x2-x1);
            continue
        end


        Ix(i,ix1)=(x2-p(i))/(x2-x1);
        Ix(i,ix2)=(p(i)-x1)/(x2-x1);
        % Linear interpolant matrix

        dIxdp(i,ix1)=-1/(x2-x1);
        dIxdp(i,ix2)=1/(x2-x1);
        % Derivative of linear interpolant matrix

     
end
    

%% Interpolate in z direction

% Assume z is constant, i.e. z(1)=z(2)=...z(end) 

    Iz=zeros(1,N(1));      % initialise z-interpolant
   

        [~,iz]=min(abs(z-zr(1))); 
        % find the index of the closest z node iz

        z1=z(iz);
        % find value of x at that node

        % find the node on the other side of p
        if z1>zr(1)            % if the known node is below zr
            z2=z1;             % rename node below and index
            iz2=iz;
            iz1=iz-1;          % define node above and index
            z1=z(iz-1); 
        elseif z1<zr(1)        % else if known node is above
            iz2=iz+1;          % define right node and index
            z2=z(iz2);
            iz1=iz;            % define above index
        elseif z1==zr(1)       % if receiver is on node (discnty problem) 
            iz1=iz;            
            Iz(1,iz1)=1;       % Interp matrix is value one at node    
            if z1==z(end)     % if zr is at last node at the bottom
                z2=z1;        % define dRdp as derivative above  
                iz2=iz1;
                iz1=iz2-1;
                z1=z(iz1);
            else              % else define as deriavtive below
                iz2=iz+1;
                z2=z(iz2);
            end
        end


        Iz(1,iz1)=(z2-zr(1))/(z2-z1);
        Iz(1,iz2)=(zr(1)-z1)/(z2-z1);
        % Linear interpolant matrix

        % Derivative of linear interpolant matrix with respect to x
        % coordinate is zero
        
        
    R=kron(Ix, Iz);                % define R as x & z interpolant
    R=full(R);                     % Size Nr x N^2
     
    dRdp=kron(dIxdp,Iz);         
    dRdp=full(dRdp);
    
end