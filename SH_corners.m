run corner_inputs.m

%% solvers ----------
% the functions "stek_helm_corners" and "stek_lap_corners" approximate the
% Steklov-Helmholtz and Steklov-Laplace eigenpairs. In case some boundary
% pieces have other boundary conditions, then the spectrum is not the pure
% Steklov spectrum. Note that you can select which eigenfunction to view in
% the "eigenfunction plots" section below. 

% outputs are as follows for the Helmholtz problem:                       %
% evs_sh = steklov_helmholtz eigenvalues                                  %
% eden_sh = steklov_helmholtz eigendensities                              %
% evecs_sh = steklov_helmholtz eigenfunctions, on the boundary            %
% A_sh = normal derivative of the single layer (matrix), this is the 
%                               adjoint of the double layer plus identity %
% B_sh = single layer (matrix)                                            %
[~,evs_sh,eden_sh,evecs_sh,A_sh,B_sh] = stek_helm_corners(curve_name,curve_params,mu,N,M,p,al,bt,len1,tol,Hom,rho_wt,eigsk,neumann_pieces,dirichlet_pieces);

% outputs are as follows for the Laplace problem:                         %
% evs_sl = steklov_laplace eigenvalues                                    %
% eden_sl = steklov_laplace eigendensities                                %
% evecs_sl = steklov_laplace eigenfunctions, on the boundary              %
% A_sl = normal derivative of the single layer (matrix), this is the 
%                               adjoint of the double layer plus identity %
% B_sl = single layer (matrix)                                            %
[len,evs_sl,eden_sl,evecs_sl,A_sl,B_sl] = stek_lap_corners(curve_name,curve_params,N,M,p,al,bt,len1,tol,Hom,rho_wt,eigsk,neumann_pieces,dirichlet_pieces);

% remove large and near 0 eigenpairs. 
% In the case of solving a "mixed boundary conditions" problem, some
% eigenvalues have really high and near 0 (or 0) magnitude. Here we remove
% them in postprocessing. 
if sum(neumann_pieces) ~= 0 || sum(dirichlet_pieces) ~= 0
inds1 = find(abs(evs_sh)> 1e4);
inds2 = find(abs(evs_sh)< 1e-3 & abs(evs_sh)> 1e-17);
inds = [inds1; inds2];
evs_sh(inds) = [];
eden_sh = eden_sh(:,setdiff(1:(2*N-length(remove)),inds));
evecs_sh = evecs_sh(:,setdiff(1:(2*N-length(remove)),inds));

inds1 = find(abs(evs_sl)> 1e4);
inds2 = find(abs(evs_sl)< 1e-3 & abs(evs_sl)> 1e-17 &abs(evs_sl)~= 0);
inds = [inds1; inds2];
evs_sl(inds) = [];
eden_sl = eden_sl(:,setdiff(1:(2*N-length(remove)),inds));
evecs_sl = evecs_sl(:,setdiff(1:(2*N-length(remove)),inds));
end
evs_sh(abs(evs_sh) > 1e4) = [];
evs_sl(abs(evs_sl) > 1e4) = [];

%% eigenvalue plots ------ 
km = min([N/2,length(evs_sh),length(evs_sl)]);
k = 1:km; % which evs you want to see ?
figure
plot(k,evs_sl(k),'--*')
hold on
plot(k,evs_sh(k),'--*')
xlabel("$k$"); ylabel("$\sigma_k$")
legend(["SL","SH"],'Interpreter','latex')

cn = strrep(curve_name,"_","\_");
title_text = join([cn,", $\mu=$ ",num2str(mu)],"");
title(title_text)

% in case of semi-circle with Steklov BC on the curve and Neumann on the
% flat part, the spectrum is the same as that of pure Steklov for the disk. 
% This is true for Dirichlet on the flat part as well.  
if sum(neumann_pieces~=0)  && curve_name == "semicircle"
hm = 100;
figure
subplot(1,2,1)
plot(evs_ac(1:hm),'--*')
hold on
plot(evs_sh(1:hm),'--*')
legend(["true","approx"])
subplot(1,2,2)
semilogy(abs(evs_ac(1:hm)-evs_sh(1:hm))./abs(evs_ac(1:hm)),'--*')
hold on
semilogy(abs(evs_ac(1:hm)-evs_sh(1:hm)),'--*')
legend(["AE","RE"])
end    

%% eigenfunction plots ------ 
% addtional input: eigen number, select integers "k" (less than N) to see
% the "k"th eigen function. Store k in variable "fixes".

% sample: this selects the smallest, first positive and 25th
% fixes = [1, find(evs_sh > 0,1), 25];  

% sample: this randomly chooses an integer between 1 and 25.
fixes = randi([1 25]);

[evec_in_sh,Xsh,Ysh] = efn_in_corners2(fixes,Co,indices,points_per_piece,dTaus,eden_sh,x,dx,points_inside,mu,M,"h");
[evec_in_sl,Xsl,Ysl] = efn_in_corners2(fixes,Co,indices,points_per_piece,dTaus,eden_sl,x,dx,points_inside,mu,M,"l");

ang = 90; cm_sh = "turbo"; cm_sl = "parula";
for i = 1:length(fixes)    
sgt_text = join([cn,", $|M|$ = ",len,", $\mu=$ ",num2str(mu)],"");

x1 = x(:,dont_remove);
efnsl = evecs_sl(:,fixes(i));
efnsh = evecs_sh(:,fixes(i));
efinsl = evec_in_sl{i};
efinsh = evec_in_sh{i};

% abs ones are finicky, vary last element for better scaling
data_scale_sl_abs = [1 1 max(abs(efinsl(efinsl~=0)))];
data_scale_sh_abs = [1 1  max(abs(efinsh(efinsh~=0)))];

data_scale_sl_real = [1 1 max(abs(real(efinsl(real(efinsl)~=0))))];
data_scale_sh_real = [1 1 max(abs(real(efinsh(real(efinsh)~=0))))];

figure()
ax1 = subplot(length(fixes),2,i);
surf(Xsl{i},Ysl{i},real(evec_in_sl{i}));
hold on
plot3(x1(1,:),x1(2,:),real(efnsl),'k','LineWidth',3)
shading interp; view(0,ang); colormap(ax1,cm_sl); %colorbar
title(join(["$\sigma^{SL}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sl(fixes(i))) ],"") )
xticks(nan); yticks(nan) 
daspect(data_scale_sl_real)

ax2 = subplot(length(fixes),2,i+length(fixes));
surf(Xsh{i},Ysh{i},real(evec_in_sh{i}));
hold on
plot3(x1(1,:),x1(2,:),real(efnsh),'k','LineWidth',3)
shading interp; view(0,ang); colormap(ax2,cm_sh); %colorbar
title(join(["$\sigma^{SH}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sh(fixes(i))) ],"") )
xticks(nan); yticks(nan) 
daspect(data_scale_sh_real)

sgtitle(join([sgt_text,", $Re(u)$"],""),'fontsize',23)
pltn = join(["figs/",curve_name,"_mu",mu,"_k",fixes(i),"_real"],"");
pltn = strrep(pltn,".","_");
% save_plot(pltn,25)

figure()
ax1 = subplot(length(fixes),2,i);
surf(Xsl{i},Ysl{i},abs(evec_in_sl{i}).^2);
hold on 
plot3(x1(1,:),x1(2,:),abs(efnsl).^2,'k','linewidth',5)
shading interp; view(0,ang); colormap(ax1,cm_sl); %colorbar
title(join(["$\sigma^{SL}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sl(fixes(i))) ],"") )
xticks(nan); yticks(nan) 
daspect(data_scale_sl_abs)

ax2 = subplot(length(fixes),2,i+length(fixes));
surf(Xsh{i},Ysh{i},abs(evec_in_sh{i}).^2);
hold on
plot3(x1(1,:),x1(2,:),abs(efnsh).^2,'k','linewidth',5)
shading interp; view(0,ang); colormap(ax2,cm_sh); %colorbar
title(join(["$\sigma^{SH}_{",num2str(fixes(i)),"}=$ ",num2str(evs_sh(fixes(i))) ],"") )
xticks(nan); yticks(nan) 
daspect(data_scale_sh_abs)
sgtitle(join([sgt_text,", $|u|^2$"],""),'fontsize',23)
pltn = strrep(pltn,"real","abssq");
% save_plot(pltn,25)
end