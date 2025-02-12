% Author: Kshitij Patil, kap15@sfu.ca %
function [A,B] = layer_pots2(mu,Tau,Tau1,T,T1,dTau,Xtau,Xt,Ytau,Yt,dXtau,dXt,dYtau,dYt,d2Xt,d2Yt)
N = size(Tau,1);
%% R, X for all pairs t, tau
X = dYt.*(Xtau-Xt)-dXt.*(Ytau-Yt);
R = sqrt((Xt-Xtau).^2+(Yt-Ytau).^2);
XR_ratio = X./R;
% ----------------------------------------
%% L, M ---------------------------------------
tic
dp1norm = sqrt(dXt.^2+dYt.^2);
dp2norm = sqrt(dXtau.^2+dYtau.^2);
dpnorm_ratio = dp2norm./dp1norm;
logsin_term = log(4*sin((T1-Tau1)/2).^2);
% logsin_term = log(4*sin((T-Tau)/2).^2);
cL = (1i*mu/2)*dpnorm_ratio;
cL1 = -mu/(2*pi)*dpnorm_ratio;
L = cL.*XR_ratio.*besselh(1,mu*R); % note that besselh(1,0) blows up so when R = 0 it's useless.

if sum(sum(isnan(L))) > N
L = L - diag(diag(L)) + diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2));
L(isnan(L)) = 0;
else
L(isnan(L)) = diag((dYt.*d2Xt-dXt.*d2Yt)./(2*pi*dp1norm.^2));
end    
L1 = cL1.*XR_ratio.*besselj(1,mu*R);
L1(isnan(L1)) = 0;

L2 = L-L1.*logsin_term;
% etl = 2*diag(log(dTau.*abs(dp2norm))).*diag(L1); 
% etl = 2*diag(log(dTau).*(dp2norm)).*diag(L1); 
etl = 2*diag(log(dTau)).*diag(L1); 
etl(isinf(etl)) = 0; etl(isnan(etl)) = 0;
L2(isnan(L2)) = diag(L)+etl;

M = (1i/2)*dp2norm.*besselh(0,mu*R);
sm = sum(sum(isnan(M)));
if sm > N
    M(isnan(M)) = 0;
end
M1 = -1/(2*pi)*dp2norm.*besselj(0,mu*R);
M2 = M-M1.*logsin_term;
mc = 1i/2-double(eulergamma)/pi;
% etm = 2*diag(log(dTau.*abs(dp2norm))).*diag(M1); 
% etm = 2*diag(log(dTau).*(dp2norm)).*diag(M1); 
etm = 2*diag(log(dTau)).*diag(M1);
etm(isinf(etm)) = 0;
% et = 0;
if sm > N
M2(isinf(M2)) = diag(dp1norm.*(mc-log(mu/2*dp1norm)/pi))+etm;
else
M2(isnan(M2)) = diag(dp1norm.*(mc-log(mu/2*dp1norm)/pi))+etm;    
end    
%-------------------------------------------

N1 = N/2;
t_scaled = (0:N-1).*T1; % the first column is 0*t, second is 1*t, so on.
cos_vecs = cos(t_scaled); % l is fixed in columns, t_j is fixed in rows.
sin_vecs = sin(t_scaled);
%L1 tp
L1_tp = L1;
% L1_tp = ftp(L1,N1,cos_vecs,sin_vecs);
% L2 tp
L2_tp = L2;
% L2_tp = ftp(L2,N1,cos_vecs,sin_vecs);
%M1 tp
M1_tp = M1;
% M1_tp = ftp(M1,N1,cos_vecs,sin_vecs);
% M2 tp
M2_tp = M2;
% M2_tp = ftp(M2,N1,cos_vecs,sin_vecs);


%% RN weights (quadrature from ColtonKress), same for any integrand (L1,M1 for us) ---------------
% diff_t_tau = T-Tau;
diff_t_tau = T1-Tau1;
RN_mat = cos(diff_t_tau);
for k1 = 2:(N1-1)
RN_mat = RN_mat+cos(k1*diff_t_tau)/k1;
end

RN_mat = -2*pi*RN_mat/N1 - pi*cos(N1*(diff_t_tau))/N1^2;
%-------------------------------------------

%% Matrices ------------------------------
A = (RN_mat.*L1_tp+pi/N1*L2_tp).*dTau; % The second term is trapezoidal rule

B = (RN_mat.*M1_tp+pi/N1*M2_tp).*dTau; % B is the effect of the S
end
