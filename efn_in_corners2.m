function [evec_in,X,Y] = efn_in_corners2(fixes,Co,indices,points_per_piece,dTaus,eden,x,dx,points_inside,mu,M,LH)
m = M;
evec_in = cell(length(fixes),1);
X = cell(length(fixes),1);
Y = cell(length(fixes),1);
shifted_inds = cumsum(points_per_piece-Co);
shifted_inds = [[1 (shifted_inds(1:end-1)+1)];shifted_inds]';
pieces = length(dTaus);

for ff = 1:length(fixes)
fix = fixes(ff);
eden_fix = eden(:,fix);


solution_fixed_point = 0;
for p1 = 1:pieces
in_ins = indices(p1,:); in_ins = (in_ins(1)+Co):in_ins(2); 
ins2 = shifted_inds(p1,:); ins2 = ins2(1):ins2(2);
N1 = points_per_piece(p1)/2; dtau = dTaus{p1}; dtau = dtau(1,(1+Co):end);
dx_in = dx(:,in_ins); 
dp2norm = sqrt(dx_in(1,:).^2+dx_in(2,:).^2); 
[Xtauin,Xin] = meshgrid(x(1,in_ins),points_inside(1,:)); 
[Ytauin,Yin] = meshgrid(x(2,in_ins),points_inside(2,:)); 
Xindiff = Xtauin-Xin; Yindiff = Ytauin-Yin;
Rin = sqrt(Xindiff.^2+Yindiff.^2);
if lower(LH) == "h"
HankRin = besselh(0,mu*Rin); % rows correspond to fixed point inside 
elseif lower(LH) == "l"
HankRin = log(Rin);
else
disp("only Steklov-Helmholtz or Steklov-Laplace supported for now")
break
end    
Recon_integrand = HankRin.*transpose(eden_fix(ins2)).*dp2norm.*dtau; % make this a tp and then sum and scale
% recint_cos1 = Recon_integrand*cos_vecs(1:end-1,1:N+1)/N;
% recint_cos1(:,[1,end]) = recint_cos1(:,[1,end])/2;
% recint_sin1 = Recon_integrand*sin_vecs(1:end-1,2:N)/N;
% recint_cos = recint_cos1*cos_vecs(:,1:N+1)';
% recint_sin = recint_sin1*sin_vecs(:,2:N)';
% recint = recint_cos+ recint_sin; recint = recint(:,1:end-1);
if lower(LH) == "h"
solution_fixed_point = solution_fixed_point+1i*pi/(4*N1)*sum(Recon_integrand,2); % trapz rule, by matlab below
elseif lower(LH) == "l"
solution_fixed_point = solution_fixed_point-1/(2*N1)*sum(Recon_integrand,2);
else
disp("only Steklov-Helmholtz or Steklov-Laplace supported for now")
break
end    

end


xlin = linspace(min(points_inside(1,:)),max(points_inside(1,:)),round(1*m));
ylin = linspace(min(points_inside(2,:)),max(points_inside(2,:)),round(1*m));

[Xin2,Yin2] = meshgrid(xlin,ylin);

evec_in{ff} = griddata(points_inside(1,:),points_inside(2,:),solution_fixed_point,Xin2,Yin2,'cubic');
% bdrypts = inpolygon(Xin2,Yin2,x(1,:),x(2,:));
% Xin2(~bdrypts) = NaN; Yin2(~bdrypts) = NaN;
X{ff} = Xin2; Y{ff} = Yin2;
end



end