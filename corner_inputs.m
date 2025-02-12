clear; close all; 
set(0,'defaultTextInterpreter','latex'); 
set(0,'defaultAxesFontSize',20)
set(0, 'DefaultLineLineWidth', 2);
curves_list = ["semicircle";"rectangle";"tear";"L";...
                "tri";"sector";"arc";"jeon2"];

%% ---------------------curve and point selection --------------------------
% ----------- All inputs are to  be given in this block -------------------%

% N = half the number of quadrature points                                 %
N = 360; % make sure 2*N is divisible by number of pieces 

% mu = real wave number, increase N on increasing mu                       %
mu = 2;sqrt(pi); 

% curves_list contains the names of curves currently implemented           %
% curve_number chooses the curve                                           % 
% 1 = semicircle, 2 = rectangle, 3 = tear, 4 = L, 5 = triangle, 6 = sector %
% additional curves can be added in the "lip_cruve.m" file.                %
curve_number = 1;
curve_name = curves_list(curve_number); curve_name = lower(curve_name);

% curve_params sets the parameters of the curves                           %
% disk_r = radius of semicircle and sector                                 %
% theta  = internal angle of the sector 
% semicircle has 1 input, disk_r and the sector has 2 [disk_r,theta]
disk_r = 1; theta = 3*pi/2;

% rectangle needs 2 inputs, [l1 l2], l1 = half length, l2 = half width     %
l1 = 1/2; l2 = 1/2;

% tear has no inputs, so leave curve_params as is                          %

% L shape needs 4 inputs, [a1 b1 a2 b2]. A corner (1,1) is fixed.          % 
% a0 = a1+a2, b0 = b1+b2 is computed internally. The corners go as         %
% (1 1) --> (1+a1 1) --> (1+a1 1+b1) --> (1+a0 1+b1) 
%                                              --> (1+a0 1+b0) --> (1 1+b0)%
a1 = 1; b1 = 1; a2 = 1; b2 = 1; 

% triangles need 3 inputs, [ta1 ta2 tb2]. 
% (1,1) is fixed. The corners go as
% (1,1) --> (ta1,1) --> (ta2,tb2)
ta1 = 2; ta2 = 3/2; tb2 = sqrt(3)/2+1;

% sample curve_params 
curve_params = disk_r; % semi-circle
% curve_params = [l1 l2]; % square of side 1
% curve_params = [a1 b1 a2 b2]; % L shape 
% curve_params = [ta1 ta2 tb2]; % triangle
% curve_params = [disk_r,theta]; % sector

% set different boundary conditions for different boundary pieces 
neumann_pieces = [2]; % neumann BC on 2nd boundary
dirichlet_pieces = []; % no dirichlet BC anywhere

% domain scaling, if perimeter is "a", make it Hom*a
Hom = 1; 

% 0 for standard length, change to any number "a" will make perimeter length "a"
len1 = 0; 

% mesh grading parameter. al (alpha) is the degree of the graded mesh.
% Higher is better, numerically tested up to 8 for some curves. 
% Too high may render the matrices highly ill-conditioned. 
al = 6; 

%inputs over --------------------------------------------------------

%% curve stuff ----------
m = 200; M = m; % points for inside the domain, required in 
% postprocessing step to reconstruct eigenfunctions

bt = 0; % keep bt 0, work in progress
p = 2^bt*al*(1+bt*log(2));

% boundary weight, can include a positive and bounded boundary weight
% function, rho_wt = @(t) sin(t)
rho_wt = 1; 
eigsk = 0; % 0 for direct solve, set positive integer "k" for first "k" eigenvalues 

tol = 1e-4; % work in progress

% cut off, keep co = 1, work in progress 
co = 0.5;
eta = 1 - 2/((1 + bt)*2 + co);
if bt > 0
Co = floor((2*N)^eta);
else
Co = 1;
end    

a = 0; b = 2*pi;

g = @(s) (s.^al).^(1./(s.^bt));
dg = @(s) al*g(s)./(s.^(bt+1)).*(1-bt*log(s));
v1 = @(s) (1/p-0.5)*(1-2*s/b).^3+(2*s/b-1)/p+0.5;
dv1 = @(s) 3*(1/p-0.5)*(-2/b)*(1-2*s/b).^2+2/(p*b);

w = @(s) b*g(v1(s))./(g(v1(s))+g(1-v1(s)));
dw = @(s) b*( dg(v1(s)).*g(1-v1(s))+g(v1(s)).*dg(1-v1(s)) ).*dv1(s)./(g(v1(s))+g(1-v1(s))).^2;

[x,dx,d2x,nx,len,dom_area,ts,ts1,Dws,indices,points_per_piece,pieces] = lip_curve(curve_name,N,p,al,bt,len1,curve_params);
indices_begin = indices(:,1); indices_end = indices(:,2);
remove0 = [1; indices_end(1:end-1)+1];
remove1 = [Co; indices_end(1:end-1)+Co];
rmv = arrayfun(@(x,y) x:y, remove0,remove1,'UniformOutput',false);
remove = [rmv{:}];

if pieces == 1
    remove = sort([remove, abs(remove-2*N)]);
end
dont_remove = setdiff(1:2*N,remove);

% dom_area = polyarea(x(1,:)',x(2,:)');
kappa = (dx(1,1:end-1).*d2x(2,1:end-1)-d2x(1,1:end-1).*dx(2,1:end-1))./(dx(1,1:end-1).^2+dx(2,1:end-1).^2).^1.5;
kappa = sum(kappa)*pi/N; kappa = kappa/(2*pi);

xi = linspace(min(x(1,:)),max(x(1,:)),m);
yi = linspace(min(x(2,:)),max(x(2,:)),m);
[X,Y] = meshgrid(xi,yi);
idx = inpolygon(X(:),Y(:),x(1,:),x(2,:)) ;
points_inside = transpose([X(idx),Y(idx)]);
%points_inside = transpose([X(:),Y(:)]);

points_inside_close = 0.99*x(:,1:end-1);
points_inside = [points_inside,points_inside_close];

box_dims = ...
    [min(x(1,:))-abs(min(x(1,:))/2) max(x(1,:))+abs(max(x(1,:))/2)...
     min(x(2,:))-abs(min(x(2,:))/2) max(x(2,:))+abs(max(x(2,:))/2)];% min x, max x,min y, max y
if box_dims(1) == 0
   box_dims(1) = -box_dims(2);
end
if box_dims(2) == 0
   box_dims(2) = -box_dims(1);
end
if box_dims(3) == 0
   box_dims(3) = -box_dims(4);
end
if box_dims(4) == 0
   box_dims(4) = -box_dims(3);
end

xi = linspace(box_dims(1),box_dims(2),round(m/2));
yi = linspace(box_dims(3),box_dims(4),round(m/2));
[X,Y] = meshgrid(xi,yi);

idx = inpolygon(X(:),Y(:),X(:),Y(:)) ;
points_inside = transpose([X(idx),Y(idx)]);

% ms = 200;
% msy = (max(x(2,1:end-1))-min(x(2,1:end-1)))/ms;
% msx = (max(x(1,1:end-1))-min(x(1,1:end-1)))/ms;
% points_inside = [min(x(1,1:end-1)):msx:max(x(1,1:end-1));min(x(2,1:end-1)):msy:max(x(2,1:end-1))];

[Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
[Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

% fig=figure(); 
hold on; 
plot(x(1,:),x(2,:),'*'); 
quiver(x(1,:),x(2,:),nx(1,:),nx(2,:),'LineWidth',2.5,'MaxHeadSize',0.5,'color',"#77AC30")
axis equal

[Taus, Ts] = cellfun(@(x) meshgrid(x), ts, 'UniformOutput', false);
[Taus1, Ts1] = cellfun(@(x) meshgrid(x), ts1, 'UniformOutput', false);
[dTaus, dTs] = cellfun(@(x) meshgrid(x), Dws, 'UniformOutput', false);

% semi-circle SN (circle S without mult.) 
if mu > 0
    evs_ac = arrayfun(@(n) (mu)*0.5*(besselj(n-1,(mu*disk_r))-besselj(n+1,(mu*disk_r)))./besselj(n,(mu*disk_r)),0:N);
    evs_ac = transpose(evs_ac');
else
    evs_ac = arrayfun(@(n) abs(mu)*0.5*(besseli(n-1,abs(mu)*len/(2*pi))+besseli(n+1,abs(mu)*len/(2*pi)))./besseli(n,abs(mu)*len/(2*pi)),0:N);
end    
evsac_mat = [evs_ac',(0:N)'];
[evs_ac,ind2] = sort(evs_ac);
evs_ac = evs_ac';
evs_ac(isinf(evs_ac)) = [];
evs_ac(isnan(evs_ac)) = [];
evs_ac(evs_ac==0) = [];
