function [A_tilde,B_tilde] = ...
    talking_pots(mu,dTau,x_in,dXtau_in,dYtau_in,x_out,dXt_out,dYt_out)

N = size(dTau,1); N1 = N/2;

% x = x_in; dXtau = dXtau_in; dYtau = dYtau_in; 
% required_points = x_out; dXt = dXt_out; dYt = dYt_out; 

x = x_in; dXtau = dXtau_in; dYtau = dYtau_in; 
dx_in = dXtau_in(1,:); dy_in = dYtau_in(1,:); 
required_points = x_out; dXt = dXt_out; dYt = dYt_out; 
dx_out = dXt(1,:); dy_out = dYt(1,:); dtau = dTau(1,:);

A_tilde = zeros(size(x,2),size(required_points,2));
B_tilde = A_tilde;
for i = 1:size(A_tilde,1)
xfix = x(1,i); yfix = x(2,i); dxfix = dx_in(:,i); dyfix = dy_in(:,i);
dp1norm = sqrt(dxfix^2+dyfix^2);
xdiff = required_points(1,:)-xfix; 
ydiff = required_points(2,:)-yfix;
Rin = sqrt(xdiff.^2+ydiff.^2); 
X = dyfix*xdiff-dxfix*ydiff;
dp2norm = sqrt(dx_out.^2+dy_out.^2); 
dpnorm_ratio = dp2norm./dp1norm; XR_ratio = X./Rin;

A2_tilde = (1i*mu/2)*dpnorm_ratio.*XR_ratio.*besselh(1,mu*Rin);
A_tilde(i,:) = pi/N1*A2_tilde.*dtau;

B2_tilde = (1i/2)*besselh(0,mu*Rin).*dp2norm; % make this a tp and then sum and scale
B_tilde(i,:) = pi/N1*B2_tilde.*dtau;
end    
% [Xtauin,Xin] = meshgrid(x(1,:),required_points(1,:)); % mesh x coord of req. points and bdry points
% [Ytauin,Yin] = meshgrid(x(2,:),required_points(2,:)); % mesh y coord of req. points and bdry points

% [Xtauin,Xin] = meshgrid(required_points(1,:),x(1,:)); % mesh x coord of req. points and bdry points
% [Ytauin,Yin] = meshgrid(required_points(2,:),x(2,:)); % mesh y coord of req. points and bdry points
% 
% Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin; % tau corresponds to the boundary under consideration
% Rin = sqrt(Xindiff.^2+Yindiff.^2); % rows correspond to required point 
% dp1norm = sqrt(dXtau.^2+dYtau.^2);
% dp2norm = sqrt(dXt.^2+dYt.^2);
% dpnorm_ratio = dp2norm./dp1norm;
% 
% X = dYt.*Xindiff-dXt.*Yindiff; % outward normal for outer boundary dotted with distance between point on both boundaries
% XR_ratio = X./Rin;
% A2_tilde = (1i*mu/2)*dpnorm_ratio.*XR_ratio.*besselh(1,mu*Rin);
% A_tilde = pi/N1*A2_tilde.*dTau;
% 
% B2_tilde = (1i/2)*besselh(0,mu*Rin).*dp2norm; % make this a tp and then sum and scale
% B_tilde = pi/N1*B2_tilde.*dTau;
end