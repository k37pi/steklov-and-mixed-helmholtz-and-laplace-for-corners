function [len,evs,eden,evecs,A,B] = stek_helm_corners(curve_name,curve_params,mu,N,M,p,al,bt,len1,tol,Hom,rho_wt,eigsk,neumann_pieces,dirichlet_pieces)

% points_per_piece = 2*N/pieces;
% indices_end = (1:pieces)*points_per_piece;
% indices_begin = (0:pieces-1)*points_per_piece+1;
% indices = [indices_begin',indices_end'];

% len1 = 0;
m = 200;
%--------------------------------------------------------

%% curve stuff ----------
% a = 0; b = 2*pi;
% dt = (b-a)/points_per_piece;
% t = a:dt:b; % 2N+1 points
% p = 5; 
% v1 = @(s) (1/p-0.5)/pi^3*(pi-s).^3+(s-pi)/(p*pi)+0.5;
% dv1 = @(s) -3*(1/p-0.5)/pi^3*(pi-s).^2+1/(p*pi);
% w = @(s) 2*pi*v1(s).^p./(v1(s).^p+v1(2*pi-s).^p);
% dw = @(s) 2*pi*(p*v1(s).^(p-1).*dv1(s).*(v1(s).^p+v1(2*pi-s).^p)-v1(s).^p.*(p*v1(s).^(p-1).*dv1(s)-p*v1(2*pi-s).^(p-1).*dv1(2*pi-s)))./(v1(s).^p+v1(2*pi-s).^p).^2;
% t1 = t; t = w(t); t2 = t;
% t((end+1)/2+1:end) = -t((end+1)/2-1:-1:1)+2*pi;

co = 0.5; % this is important, increase with N
eta = 1 - 2/((1 + bt)*2 + co);
if bt > 0
Co = floor((2*N)^eta);
else
Co = 1;
end

% tp = 2*N;
% a = 0; b = 2*pi;
g = @(s) (s.^al).^(1./(s.^bt));
dg = @(s) al*g(s)./(s.^(bt+1)).*(1-bt*log(s));
v1 = @(s) (1/p-0.5)*(1-2*s/b).^3+(2*s/b-1)/p+0.5;
dv1 = @(s) 3*(1/p-0.5)*(-2/b)*(1-2*s/b).^2+2/(p*b);

% if bt > 0 
w = @(s) b*g(v1(s))./(g(v1(s))+g(1-v1(s)));
dw = @(s) b*( dg(v1(s)).*g(1-v1(s))+g(v1(s)).*dg(1-v1(s)) ).*dv1(s)./(g(v1(s))+g(1-v1(s))).^2;
% else
% w = @(s) 2*pi*v1(s).^p./(v1(s).^p+v1(2*pi-s).^p);
% dw = @(s) 2*pi*(p*v1(s).^(p-1).*dv1(s).*(v1(s).^p+v1(2*pi-s).^p)-v1(s).^p.*(p*v1(s).^(p-1).*dv1(s)-p*v1(2*pi-s).^(p-1).*dv1(2*pi-s)))./(v1(s).^p+v1(2*pi-s).^p).^2;
% end

[x,dx,d2x,~,len,~,ts,ts1,Dws,indices,points_per_piece,pieces] = lip_curve(curve_name,N,p,al,bt,len1,curve_params);
points_per_piece = round(points_per_piece);
indices_begin = indices(:,1); indices_end = indices(:,2);

% points_per_piece = indices(:,2)-indices(:,1)+1;
% remove = [1; indices_end(1:end-1)+1];
% remove = [1:Co; (indices_end(1:end-1)+1):(indices_end(1:end-1)+Co)];
remove0 = [1; indices_end(1:end-1)+1];
remove1 = [Co; indices_end(1:end-1)+Co];
rmv = arrayfun(@(x,y) x:y, remove0,remove1,'UniformOutput',false);
remove = [rmv{:}];

if pieces == 1
    remove = sort([remove, abs(remove-2*N)]);
end

% [x,dx,d2x,nx,len] = lip_curve(curve_name,t,len1,curve_params);
% [x,dx,d2x,~,len,~,ts,ts1,Dws,indices,points_per_piece,pieces] = lip_curve(curve_name,N,p,al,bt,len1,curve_params);
% indices_begin = indices(:,1); indices_end = indices(:,2);
% % points_per_piece = indices(:,2)-indices(:,1)+1;
% remove = [1; indices_end(1:end-1)+1];
dont_remove = setdiff(1:2*N,remove);

xi = linspace(min(x(1,:)),max(x(1,:)),m);
yi = linspace(min(x(2,:)),max(x(2,:)),m);
[X,Y] = meshgrid(xi,yi);
idx = inpolygon(X(:),Y(:),x(1,:),x(2,:)) ;
% points_inside = transpose([X(idx),Y(idx)]);
%points_inside = transpose([X(:),Y(:)]);

% points_inside_close = 0.99*x(:,1:end-1);
% points_inside = [points_inside,points_inside_close];

% ms = 200;
% msy = (max(x(2,1:end-1))-min(x(2,1:end-1)))/ms;
% msx = (max(x(1,1:end-1))-min(x(1,1:end-1)))/ms;
% points_inside = [min(x(1,1:end-1)):msx:max(x(1,1:end-1));min(x(2,1:end-1)):msy:max(x(2,1:end-1))];

% [Xtauin,Xin] = meshgrid(x(1,1:end-1),points_inside(1,:)); % mesh x coord of in points and bdry points
% [Ytauin,Yin] = meshgrid(x(2,1:end-1),points_inside(2,:)); % mesh y coord of in points and bdry points

[Taus, Ts] = cellfun(@(x) meshgrid(x), ts, 'UniformOutput', false);
[Taus1, Ts1] = cellfun(@(x) meshgrid(x), ts1, 'UniformOutput', false);
[dTaus, ~] = cellfun(@(x) meshgrid(x), Dws, 'UniformOutput', false);
% [Tau,T] = meshgrid(t);
% [Tau1,T1] = meshgrid(t1);
% [dTau,dT] = meshgrid(Dw);
% [Xtau,Xt] = meshgrid(x(1,:));
% [Ytau,Yt] = meshgrid(x(2,:));
% [dXtau,dXt] = meshgrid(dx(1,:));
% [dYtau,dYt] = meshgrid(dx(2,:)); 
% [d2Xtau,d2Xt] = meshgrid(d2x(1,:));
% [d2Ytau,d2Yt] = meshgrid(d2x(2,:));

%% lonely and talking pieces -------
% As = cell(pieces,1); Bs = cell(pieces,1);
A_tildes = cell(pieces,(pieces)); B_tildes = cell(pieces,(pieces));
% for pp = 1:pieces
% inds = indices_begin(pp):indices_end(pp);    
% [Tau,T] = meshgrid(t);
% [dTau,dT] = meshgrid(dw(t));
% [Xtau,Xt] = meshgrid(x(1,inds));
% [Ytau,Yt] = meshgrid(x(2,inds));
% [dXtau,dXt] = meshgrid(dx(1,inds));
% [dYtau,dYt] = meshgrid(dx(2,inds)); 
% [d2Xtau,d2Xt] = meshgrid(d2x(1,inds));
% [d2Ytau,d2Yt] = meshgrid(d2x(2,inds));    
% [As{pp},Bs{pp}] =...
%     layer_pots2(mu,Tau,T,dTau,Xtau,Xt,Ytau,Yt,dXtau,dXt,dYtau,dYt,d2Xt,d2Yt);
% end
%----------------------------------------

for pp = 1:pieces
alls = 1:pieces; alls(pp) = [];    
inds = indices_begin(pp):indices_end(pp);    
[Xtau,Xt] = meshgrid(x(1,inds));
[Ytau,Yt] = meshgrid(x(2,inds));
[dXtau,dXt] = meshgrid(dx(1,inds));
[dYtau,dYt] = meshgrid(dx(2,inds)); 
[~,d2Xt] = meshgrid(d2x(1,inds));
[~,d2Yt] = meshgrid(d2x(2,inds));    

% lonely 
[A_tildes{pp,pp},B_tildes{pp,pp}] =...
    layer_pots2(mu,Taus{pp},Taus1{pp},Ts{pp},Ts1{pp},dTaus{pp},Xtau,Xt,Ytau,Yt,dXtau,dXt,dYtau,dYt,d2Xt,d2Yt);

for pp2 = alls
inds_out = indices_begin(pp2):indices_end(pp2);    
% [Xtau_out,Xt_out] = meshgrid(x(1,inds_out));
% [Ytau_out,Yt_out] = meshgrid(x(2,inds_out));
[~,dXt_out] = meshgrid(dx(1,inds_out));
[~,dYt_out] = meshgrid(dx(2,inds_out)); 
% [d2Xtau_out,d2Xt_out] = meshgrid(d2x(1,inds_out));
% [d2Ytau_out,d2Yt_out] = meshgrid(d2x(2,inds_out));

x_in = x(:,inds); x_out = x(:,inds_out);

% talking 
[A_tildes{pp,pp2},B_tildes{pp,pp2}] = ...
    talking_pots(mu,dTaus{pp2},x_in,dXtau,dYtau,x_out,dXt_out,dYt_out); 
% x_in is the fixed point on the boundary, x_out is the boundary being integrated
% on
end
end

% I = eye(size(Tau));

%% matrix assembly line ----------
A = zeros(2*N,2*N); B = A;
for i = 1:pieces
row_inds = indices(i,1):indices(i,2);
for j = 1:pieces    
col_inds = indices(j,1):indices(j,2);    
% A(row_inds,col_inds) = A_tildes{i,j};
% if i == j
%     A(row_inds,col_inds) = A(row_inds,col_inds) + eye(points_per_piece(i));
% end    
% % B(row_inds,col_inds) = B_tildes{i,j};
% 
% if sum(neumann_pieces==i) == 0 
% B(row_inds,col_inds) = B_tildes{i,j};
% end    
% A = [As{1}+I A_tildes{1,2}; A_tildes{2,1} As{2}+I];
% % B = [Bs{1} B_tildes{1,2}; 0*Bs{2} 0*B_tildes{2,1}];
% B = [Bs{1} B_tildes{1,2}; Bs{2} B_tildes{2,1}];
if sum(neumann_pieces==i) == 0 
B(row_inds,col_inds) = B_tildes{i,j};
end    
if sum(dirichlet_pieces==i) == 0 
    A(row_inds,col_inds) = A_tildes{i,j};
if i == j
    A(row_inds,col_inds) = A(row_inds,col_inds) + eye(points_per_piece(i));
end    
end
end
end

Rho_wt = 1;
B_wt = Rho_wt.*B;

%% Steklov eigenstuff ------------
Mu = mu;
[u,s,v] = svd(B_wt);
tol = 1e-6;
ss = diag(s);
%rank_def = sum(ss < 5*10^(tol));
% orders = floor(log10(ss));
% diffs = abs(ss-tol);
%orders(diffs < 1e-3) = orders(diffs < 1e-3)-1; 
% rank_def = sum(orders < log10(tol));
A = A(dont_remove,dont_remove);
B_wt = B_wt(dont_remove,dont_remove);
B = B(dont_remove,dont_remove);
rank_def1 = sum(abs(ss-tol) < tol);
rank_def = size(B,1)-rank(B);

if 1==2%rank_def ~= 0 && sum(neumann_pieces)==0
    %[u,s,v] = svd(B_wt);
    A_svd = u(:,1:end-rank_def)'*(A)*v(:,1:end-rank_def); 
    B_wt1 = s(1:end-rank_def,1:end-rank_def);
    
    %B_wt1 = u(:,1:end-rank_def)'*B_wt*v(:,1:end-rank_def);    
    % A_svd = u'*(A+I)*v; 
    % B_wt1 = u'*B_wt*v;
    [eden1,evs1] = eig(A_svd,B_wt1);
    %[eden1,evs0] = eig(A+I,B_wt);
    evs = real(diag(evs1));    
    [evs, ind] = sort(evs);%,"descend");
    [~,blah,ev_reps] = uniquetol(real(evs)); 
    ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals
    eden1 = v(:,1:end-rank_def)*eden1;
    evecs = B*eden1/2; %disk_r
    % normder_evecs = (A)*eden1/2;
    eden = eden1(:,ind); % eigen densities
    evecs = evecs(:,ind);
    % normder_evecs = normder_evecs(:,ind);    
    A = A_svd; B = B_wt1; %eden = eden1;
else
    % I = diag(ones(2*N,1));
    % A = (A+A.')/2;B_wt = (B_wt+B_wt.')/2; B = (B+B.')/2;
    % A = (A+A')/2;B_wt = (B_wt+B_wt')/2; B = (B+B')/2;
    [eden1,evs1] = eig(A,B_wt); 
    evs0 = diag(evs1);
    evs = real(evs0);    
    [evs, ind] = sort(evs);%,"descend");
    [~,blah,ev_reps] = uniquetol(real(evs)); 
    ev_reps = blah(accumarray(ev_reps,1).'==1);% indices of non repeated evals
    evecs = B*eden1/2; %disk_r
    % normder_evecs = (A)*eden1/2;
    eden = eden1(:,ind); % eigen densities
    evecs = evecs(:,ind);
    % normder_evecs = normder_evecs(:,ind);    
end    

end