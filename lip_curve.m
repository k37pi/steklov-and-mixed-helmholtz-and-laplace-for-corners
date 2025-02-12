% This function takes the curve name and returns
% the parameterization x, its derivatie dx and second derivative
% d2x
% Inputs : curve name, point set = t, normalize length = len1, curve params if any
% Outputs : x, dx, d2x at given point set
function [x,dx,d2x,nx,len,dom_area,ts,ts1,Dws,indices,points_per_piece,pieces] = lip_curve(name,N,p,al,bt,len1,varargin)%curve(name,t,len1,a,varargin)
    % v1 = @(s) (1/p-0.5)/pi^3*(pi-s).^3+(s-pi)/(p*pi)+0.5;
    % dv1 = @(s) -3*(1/p-0.5)/pi^3*(pi-s).^2+1/(p*pi);
    % w = @(s) 2*pi*v1(s).^p./(v1(s).^p+v1(2*pi-s).^p);
    % dw = @(s) 2*pi*(p*v1(s).^(p-1).*dv1(s).*(v1(s).^p+v1(2*pi-s).^p)-v1(s).^p.*(p*v1(s).^(p-1).*dv1(s)-p*v1(2*pi-s).^(p-1).*dv1(2*pi-s)))./(v1(s).^p+v1(2*pi-s).^p).^2;

    % al = 6; bt = 0;
    % p = 2^bt*(al+bt*log(2));
    tp = 2*N;
    a = 0; b = 2*pi;

    g = @(s) ((s).^al).^(1./((s).^bt));
    dg = @(s) al*g(s)./(s.^(bt+1)).*(1-bt*log(s));
    v1 = @(s) (1/p-0.5)*(1-2*s/b).^3+(2*s/b-1)/p+0.5;
    dv1 = @(s) 3*(1/p-0.5)*(-2/b)*(1-2*s/b).^2+2/(p*b);
    
    w = @(s) b*g(v1(s))./(g(v1(s))+g(1-v1(s)));
    dw = @(s) b*( dg(v1(s)).*g(1-v1(s))+g(v1(s)).*dg(1-v1(s)) ).*dv1(s)./(g(v1(s))+g(1-v1(s))).^2;
    
    % g = @(s) (s.^al).^(1./(s.^bt));
    % dg = @(s) al*g(s)./(s.^(bt+1)).*(1-bt*log(s));
    % 
    % v1 = @(s) (1/p-0.5)*(1-2*s/b).^3+(2*s/b-1)/p+0.5;
    % dv1 = @(s) 3*(1/p-0.5)*(-2/b)*(1-2*s/b).^2+2/(p*b);
    % 
    % if bt > 0 
    % w = @(s) b*g(v1(s))./(g(v1(s))+g(1-v1(s)));
    % dw = @(s) b*( dg(v1(s)).*g(1-v1(s))+g(v1(s)).*dg(1-v1(s)) ).*dv1(s)./(g(v1(s))+g(1-v1(s))).^2;
    % else
    % w = @(s) 2*pi*v1(s).^p./(v1(s).^p+v1(2*pi-s).^p);
    % dw = @(s) 2*pi*(p*v1(s).^(p-1).*dv1(s).*(v1(s).^p+v1(2*pi-s).^p)-v1(s).^p.*(p*v1(s).^(p-1).*dv1(s)-p*v1(2*pi-s).^(p-1).*dv1(2*pi-s)))./(v1(s).^p+v1(2*pi-s).^p).^2;
    % end
    switch lower(name)
        case 'semicircle'
            % total_bdry = pi+2;
            % p1 = pi/total_bdry ~ 0.611 ~ 0.6
            p1N = 0.6*tp; p2N = tp-p1N;
            % p1N = 0.5*tp; p2N = tp-p1N;
            pieces = 2; %ts = cell(1,pieces);
            
            points_per_piece = [p1N p2N];
            % indices_end = zeros(pieces,2); 
            indices_begin = [1 points_per_piece(1:end-1)+1];
            indices_end = cumsum(points_per_piece);
            % indices_end = (1:pieces)*points_per_piece;
            % indices_begin = (0:pieces-1)*points_per_piece+1;
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 

            r = varargin{1}; t = ts{1};
            xs = r*[cos(t/2);sin(t/2)];
            xs(:,((end+1)/2+1):end) = [-xs(1,((end+1)/2-1):-1:1);xs(2,((end+1)/2-1):-1:1) ];
            dxs = r/2*[-sin(t/2);cos(t/2)];
            dxs(:,((end+1)/2+1):end) = [dxs(1,((end+1)/2-1):-1:1);-dxs(2,((end+1)/2-1):-1:1)];
            d2xs = r/4*[-cos(t/2);-sin(t/2)];
            d2xs(:,((end+1)/2+1):end) = [-d2xs(1,((end+1)/2-1):-1:1);d2xs(2,((end+1)/2-1):-1:1)];
            xs = xs(:,1:end-1); dxs = dxs(:,1:end-1); d2xs = d2xs(:,1:end-1);
            
            t = ts{2};
            xl = r*[t/pi-1;0*t];
            xl(:,((end+1)/2+1):end) = [-xl(1,((end+1)/2-1):-1:1);xl(2,((end+1)/2-1):-1:1) ];
            dxl = r*[1/pi*(t-t+1);0*t];
            d2xl = [0*t;0*t];
            xl = xl(:,1:end-1); dxl = dxl(:,1:end-1); d2xl = d2xl(:,1:end-1);
            x = [xs,xl]; 
            dx = [dxs,dxl];
            d2x = [d2xs,d2xl];
            % x(1,t >= pi) = r*(2*t(t >= pi)/pi-3);
            % x(2,t >= pi) = zeros(1,length(t(t >= pi)));
            % dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            % dx(1,t >= pi) = 2*r*ones(1,length(t(t >= pi)))/pi;
            % dx(2,t >= pi) = zeros(1,length(t(t >= pi)));
            % d2x = -x;
            % d2x(:,t >= pi) = zeros(2,length(t(t >= pi)));
            len = (pi+2)*r;
            % SHOULD MATCH len
            % 2*pi/length(t)*sum(sqrt(dx(1,:).^2+dx(2,:).^2)); 
            % 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end/2).^2+dx(2,1:end/2).^2))+2*pi/(length(t)-1)*sum(sqrt(dx(1,end/2+1:end).^2+dx(2,end/2+1:end).^2))
        case 'arc'
            % total_bdry = pi+2;
            % p1 = pi/total_bdry ~ 0.611 ~ 0.6
            % p1N = 0.6*tp; p2N = tp-p1N;
            % p1N = 0.5*tp; p2N = tp-p1N;
            pieces = 1; %ts = cell(1,pieces);
            
            points_per_piece = tp;
            % indices_end = zeros(pieces,2); 
            indices_begin = [1 points_per_piece(1:end-1)+1];
            indices_end = cumsum(points_per_piece);
            % indices_end = (1:pieces)*points_per_piece;
            % indices_begin = (0:pieces-1)*points_per_piece+1;
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 

            r = varargin{1}; t = ts{1};
            xs = r*[cos(t/2);sin(t/2)];
            xs(:,((end+1)/2+1):end) = [-xs(1,((end+1)/2-1):-1:1);xs(2,((end+1)/2-1):-1:1) ];
            dxs = r/2*[-sin(t/2);cos(t/2)];
            dxs(:,((end+1)/2+1):end) = [dxs(1,((end+1)/2-1):-1:1);-dxs(2,((end+1)/2-1):-1:1)];
            d2xs = r/4*[-cos(t/2);-sin(t/2)];
            d2xs(:,((end+1)/2+1):end) = [-d2xs(1,((end+1)/2-1):-1:1);d2xs(2,((end+1)/2-1):-1:1)];
            xs = xs(:,1:end-1); dxs = dxs(:,1:end-1); d2xs = d2xs(:,1:end-1);
            
            % t = ts{2};
            % xl = r*[t/pi-1;0*t];
            % xl(:,((end+1)/2+1):end) = [-xl(1,((end+1)/2-1):-1:1);xl(2,((end+1)/2-1):-1:1) ];
            % dxl = r*[1/pi*(t-t+1);0*t];
            % d2xl = [0*t;0*t];
            % xl = xl(:,1:end-1); dxl = dxl(:,1:end-1); d2xl = d2xl(:,1:end-1);
            x = xs; 
            dx = dxs;
            d2x = d2xs;
            % x(1,t >= pi) = r*(2*t(t >= pi)/pi-3);
            % x(2,t >= pi) = zeros(1,length(t(t >= pi)));
            % dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            % dx(1,t >= pi) = 2*r*ones(1,length(t(t >= pi)))/pi;
            % dx(2,t >= pi) = zeros(1,length(t(t >= pi)));
            % d2x = -x;
            % d2x(:,t >= pi) = zeros(2,length(t(t >= pi)));
            len = pi*r;
            % SHOULD MATCH len
            % 2*pi/length(t)*sum(sqrt(dx(1,:).^2+dx(2,:).^2)); 
            % 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end/2).^2+dx(2,1:end/2).^2))+2*pi/(length(t)-1)*sum(sqrt(dx(1,end/2+1:end).^2+dx(2,end/2+1:end).^2))
        case 'sector'
            % total_bdry = pi+2;
            % p1 = pi/total_bdry ~ 0.611 ~ 0.6
            ab = varargin{1};
            r = ab(1);theta = ab(2);
            ratios = [1 theta 1]/(2+theta);
            ratios = round(ratios*10); 
            if sum(ratios) > 10
            ratios(ratios==(max(ratios))) = max(ratios)-1;
            end
            if sum(ratios) < 10
            ratios([1 3]) = (10-max(ratios))/2+zeros(1,2);
            end
            ratios = ratios/10;
            % p1N = 0.6*tp; p2N = tp-p1N;
            % p1N = 0.5*tp; p2N = tp-p1N;
            pieces = 3; %ts = cell(1,pieces);
            xi = (pi-theta)/2;
            
            points_per_piece = ratios*tp;
            indices_end = cumsum(points_per_piece);
            indices_begin = [1 indices_end(1:end-1)+1];
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 
            
            t = ts{1};
            spx = 1; spy1 = 1; epx = 1+r*cos(xi); epy = 1+r*sin(xi);
            [xl1, dxl1, d2xl1] = slantline(t,spx,epx,spy1,epy); 
            
            epx = 1; epy = 1; spx = 1+r*cos(xi+theta); spy1 = 1+r*sin(xi+theta);
            t = ts{3};
            [xl2, dxl2, d2xl2] = slantline(t,spx,epx,spy1,epy); 

            t = ts{2};
            xs = [1+r*cos(xi+t*theta/(2*pi));1+r*sin(xi+t*theta/(2*pi))];
            mp = xs(1,(end+1)/2);
            xs(1,((end+1)/2+1):end) = 2*mp-xs(1,((end+1)/2-1):-1:1);
            xs(2,((end+1)/2+1):end) = xs(2,((end+1)/2-1):-1:1);

            dxs = r*theta/(2*pi)*[-sin(xi+t*theta/(2*pi));cos(xi+t*theta/(2*pi))];
            % dxs(:,((end+1)/2+1):end) = [dxs(1,((end+1)/2-1):-1:1);-dxs(2,((end+1)/2-1):-1:1)];
            % mp = dxs(1,(end+1)/2);
            dxs(1,((end+1)/2+1):end) = dxs(1,((end+1)/2-1):-1:1);
            mp = dxs(2,(end+1)/2);
            dxs(2,((end+1)/2+1):end) = 2*mp-dxs(2,((end+1)/2-1):-1:1);

            d2xs = r*theta^2/(2*pi)^2*[-cos(xi+t*theta/(2*pi));-sin(xi+t*theta/(2*pi))];
            % d2xs(:,((end+1)/2+1):end) = [-d2xs(1,((end+1)/2-1):-1:1);d2xs(2,((end+1)/2-1):-1:1)];
            mp = d2xs(1,(end+1)/2);
            d2xs(1,((end+1)/2+1):end) = 2*mp-d2xs(1,((end+1)/2-1):-1:1);
            % mp = d2xs(2,(end+1)/2);
            d2xs(2,((end+1)/2+1):end) = d2xs(2,((end+1)/2-1):-1:1);
            
            xs = xs(:,1:end-1); dxs = dxs(:,1:end-1); d2xs = d2xs(:,1:end-1);
            
            x = [xl1,xs,xl2]; 
            dx = [dxl1,dxs,dxl2];
            d2x = [d2xl1,d2xs,d2xl2];
            % x(1,t >= pi) = r*(2*t(t >= pi)/pi-3);
            % x(2,t >= pi) = zeros(1,length(t(t >= pi)));
            % dx = r*[-sin(t);cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            % dx(1,t >= pi) = 2*r*ones(1,length(t(t >= pi)))/pi;
            % dx(2,t >= pi) = zeros(1,length(t(t >= pi)));
            % d2x = -x;
            % d2x(:,t >= pi) = zeros(2,length(t(t >= pi)));
            len = (theta+2)*r;
            % SHOULD MATCH len
            % 2*pi/length(t)*sum(sqrt(dx(1,:).^2+dx(2,:).^2)); 
            % 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end/2).^2+dx(2,1:end/2).^2))+2*pi/(length(t)-1)*sum(sqrt(dx(1,end/2+1:end).^2+dx(2,end/2+1:end).^2))
        case 'rectangle'
            ab = varargin{1}; 
            a1 = ab(1); b1 = ab(2);
            len = 4*(a1+b1); pieces = 4;
            ratios = 2*[b1 a1 b1 a1]/len;
            % ratios = [1 1 1 1]/4;
            if abs(sum(ratios)-1) > 1e-15
            % ratios(ratios==min(ratios)) = ratios(ratios==min(ratios))+1-sum(ratios); 
            disp("non proper distribution of points")
            end
            points_per_piece = ratios*tp;
            indices_end = cumsum(points_per_piece);
            indices_begin = [1 indices_end(1:end-1)+1];
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 

            t = ts{1};
            x1 = [a1+0*t;b1*(t/pi-1)];
            x1(2,((end+1)/2+1):end) = -x1(2,((end+1)/2-1):-1:1);
            dx1 = [0*t;b1/pi+0*t];
            d2x1 = [0*t;0*t];
            x1 = x1(:,1:end-1); dx1 = dx1(:,1:end-1); d2x1 = d2x1(:,1:end-1);
            
            t = ts{3};
            x3 = [-a1+0*t;-b1*(t/pi-1)];
            x3(2,((end+1)/2+1):end) = -x3(2,((end+1)/2-1):-1:1);
            dx3 = [0*t;-b1/pi+0*t];
            d2x3 = [0*t;0*t];
            x3 = x3(:,1:end-1); dx3 = dx3(:,1:end-1); d2x3 = d2x3(:,1:end-1);
            
            t = ts{2};
            x2 = [-a1*(t/pi-1);b1+0*t];
            x2(1,((end+1)/2+1):end) = -x2(1,((end+1)/2-1):-1:1);
            dx2 = [-a1/pi+0*t;0*t];
            d2x2 = [0*t;0*t];
            x2= x2(:,1:end-1); dx2 = dx2(:,1:end-1); d2x2 = d2x2(:,1:end-1);
            
            t = ts{4};
            x4 = [a1*(t/pi-1);-b1+0*t];
            x4(1,((end+1)/2+1):end) = -x4(1,((end+1)/2-1):-1:1);
            dx4 = [a1/pi+0*t;0*t];
            d2x4 = [0*t;0*t];
            x4 = x4(:,1:end-1); dx4 = dx4(:,1:end-1); d2x4 = d2x4(:,1:end-1);

            x = [x1,x2,x3,x4]; 
            dx = [dx1,dx2,dx3,dx4];
            d2x = [d2x1,d2x2,d2x3,d2x4];
            
            % SHOULD MATCH len
            % 2*pi/length(t)*sum(sqrt(dx(1,:).^2+dx(2,:).^2)); 
            % 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end/2).^2+dx(2,1:end/2).^2))+2*pi/(length(t)-1)*sum(sqrt(dx(1,end/2+1:end).^2+dx(2,end/2+1:end).^2))
        case 'tri'
            ab = varargin{1}; 
            a1 = ab(1); a2 = ab(2); b2 = ab(3);
            pieces = 3;
            % first point (1,1)
            % second point (a1,1),eg (2,1)
            % third point (a2,b2), eg (1.5,sqrt(3)/2+1)
            gl1 = a1-1; gl2 = sqrt((a2-a1)^2+(b2-1)^2);
            gl3 = sqrt((a2-1)^2+(b2-1)^2);
            len = gl1+gl2+gl3;
            ratios = [gl1 gl2 gl3]/len;

            if ratios(1)~= ratios(2)
            ratios = round(ratios*10); 
            if sum(ratios) > 10
            ratios(ratios==(max(ratios))) = max(ratios)-1;
            end
            ratios = ratios/10;
            end
            % ratios = ones(1,pieces)/pieces;
            % ratios = [1 1 1 1]/4;
            if abs(sum(ratios)-1) > 1e-15
            % ratios(ratios==min(ratios)) = ratios(ratios==min(ratios))+1-sum(ratios); 
            disp("non proper distribution of points")
            end
            points_per_piece = ratios*tp;
            indices_end = cumsum(points_per_piece);
            indices_begin = [1 indices_end(1:end-1)+1];
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 
    
            t = ts{1};
            x1 = [1+(a1-1)/(2*pi)*t;1+0*t];
            mp = x1(1,(end+1)/2);
            x1(1,((end+1)/2+1):end) = 2*mp-x1(1,((end+1)/2-1):-1:1);
            dx1 = [(a1-1)/(2*pi)+0*t;0*t];
            d2x1 = [0*t;0*t];
            x1 = x1(:,1:end-1); dx1 = dx1(:,1:end-1); d2x1 = d2x1(:,1:end-1);
            
            t = ts{2};
            x2 = [a1+(a2-a1)*t/(2*pi);1+(b2-1)*t/(2*pi)];
            mp = x2(1,(end+1)/2);
            x2(1,((end+1)/2+1):end) = 2*mp-x2(1,((end+1)/2-1):-1:1);
            mp = x2(2,(end+1)/2);
            x2(2,((end+1)/2+1):end) = 2*mp-x2(2,((end+1)/2-1):-1:1);
            dx2 = [(a2-a1)+0*t;(b2-1)+0*t]/( 2*pi);
            d2x2 = [0*t;0*t];
            x2= x2(:,1:end-1); dx2 = dx2(:,1:end-1); d2x2 = d2x2(:,1:end-1);
            
            t = ts{3};
            x3 = [a2+(1-a2)*t/(2*pi);b2+(1-b2)*t/(2*pi)];
            mp = x3(1,(end+1)/2);
            x3(1,((end+1)/2+1):end) = 2*mp-x3(1,((end+1)/2-1):-1:1);
            mp = x3(2,(end+1)/2);
            x3(2,((end+1)/2+1):end) = 2*mp-x3(2,((end+1)/2-1):-1:1);

            dx3 = [(1-a2)+0*t;(1-b2)+0*t]/( 2*pi);
            d2x3 = [0*t;0*t];
            x3 = x3(:,1:end-1); dx3 = dx3(:,1:end-1); d2x3 = d2x3(:,1:end-1);

            x = [x1,x2,x3]; 
            dx = [dx1,dx2,dx3];
            d2x = [d2x1,d2x2,d2x3];
            
            % SHOULD MATCH len
            % 2*pi/length(t)*sum(sqrt(dx(1,:).^2+dx(2,:).^2)); 
            % 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end/2).^2+dx(2,1:end/2).^2))+2*pi/(length(t)-1)*sum(sqrt(dx(1,end/2+1:end).^2+dx(2,end/2+1:end).^2))

        case 'l'
            ab = varargin{1}; 
            a1 = ab(1); b1 = ab(2); 
            a2 = ab(3); b2 = ab(4);
            a0 = a1+a2; b0 = b1+b2;
            len = 2*(a0+b0); pieces = 6;

            % a1 = 1.1; b1 = 1; a2 = 0.9; b2 = 1; 
            % a0 = a1+a2; b0 = b1+b2;
            pts = [1 1; 1+a1 1; 1+a1 1+b1; 1+a0 1+b1;1+a0 1+b0; 1 1+b0];
            pts = [pts; pts(1,:)];
            % total_bdry = 2*(a0+b0);
            % p1 = pi/total_bdry ~ 0.611 ~ 0.61
            % p1N = 0.61*tp; p2N = tp-p1N;
            % pieces = 2; %ts = cell(1,pieces);
            % 
            % ratios = (round([a1 b1 a2 b2 a0 b0]/len,2));
            ratios = [a1 b1 a2 b2 a0 b0]/len;
            if sum(ratios)~=1
            % ratios(ratios==min(ratios)) = ratios(ratios==min(ratios))+1-sum(ratios); 
            disp("non proper distribution of points")
            end
            points_per_piece = ratios*tp;
            % indices_end = zeros(pieces,2); 
            % indices_begin = [1 points_per_piece(1:end-1)+1];
            indices_end = cumsum(points_per_piece);
            indices_begin = [1 indices_end(1:end-1)+1];
            % indices_end = (1:pieces)*points_per_piece;
            % indices_begin = (0:pieces-1)*points_per_piece+1;
            indices = [indices_begin',indices_end'];

            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 

            for l = 1:6
            sp = pts(l,:); ep = pts(l+1,:);     
            t = ts{l};
            if rem(l,2) == 1
            [xh,mh] = horverline(t,sp(1),ep(1));
            x1 = [xh;ep(2)+0*t];
            dx1 = [mh+0*t;0*t];
            else
            [yv,mv] = horverline(t,sp(2),ep(2));
            x1 = [ep(1)+0*t;yv];
            dx1 = [0*t;mv+0*t];
            end
            d2x1 = [0*t;0*t];
            x1 = x1(:,1:end-1); dx1 = dx1(:,1:end-1); d2x1 = d2x1(:,1:end-1);
            
            if l == 1
            x = x1; dx = dx1; d2x = d2x1;
            else
            x = [x,x1]; dx = [dx,dx1]; d2x = [d2x,d2x1];
            end

            end

            % SHOULD MATCH len
            % 2*pi/length(t)*sum(sqrt(dx(1,:).^2+dx(2,:).^2)); 
            % 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end/2).^2+dx(2,1:end/2).^2))+2*pi/(length(t)-1)*sum(sqrt(dx(1,end/2+1:end).^2+dx(2,end/2+1:end).^2))    
        
        case 'tear'
            pieces = 1; %ts = cell(1,pieces);
            
            points_per_piece = tp;
            % indices_end = zeros(pieces,2); 
            indices_begin = [1 points_per_piece(1:end-1)+1];
            indices_end = cumsum(points_per_piece);
            % indices_end = (1:pieces)*points_per_piece;
            % indices_begin = (0:pieces-1)*points_per_piece+1;
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 
            
            t = ts{1};
            x = [2*sin(t/2);-sin(t)];
            dx = [cos(t/2);-cos(t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = [-1/2*sin(t/2);sin(t)];
            x = x(:,1:end-1); dx = dx(:,1:end-1); d2x = d2x(:,1:end-1);
            % len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));
        
        case 'jeon2'
            pieces = 1; %ts = cell(1,pieces);
            
            points_per_piece = tp;
            % indices_end = zeros(pieces,2); 
            indices_begin = [1 points_per_piece(1:end-1)+1];
            indices_end = cumsum(points_per_piece);
            % indices_end = (1:pieces)*points_per_piece;
            % indices_begin = (0:pieces-1)*points_per_piece+1;
            indices = [indices_begin',indices_end'];
            
            dts = (b-a)./(points_per_piece); %cutoff = 0.1; 
            ts = arrayfun(@(x) a:x:b,dts,'UniformOutput',false); % 2N+1 points
            ts1 = ts; ts = cellfun(@(x) w(x),ts1,'UniformOutput',false);
            Dws = cellfun(@(x)dw(x),ts1,'UniformOutput',false);%dw(t1); 
            
            a1 = 3/2; a2 = 1;
            t = ts{1};
            x = [-2/3*sin(a1*t);-sin(a2*t)];
            dx = [-2/3*a1*cos(a1*t);-a2*cos(a2*t)]; %nx = [dx(2,:);-dx(1,:)];
            d2x = [2/3*a1^2*sin(a1*t);a2^2*sin(a2*t)];
            x = x(:,1:end-1); dx = dx(:,1:end-1); d2x = d2x(:,1:end-1);
            % len = 2*pi/(length(t)-1)*sum(sqrt(dx(1,1:end-1).^2+dx(2,1:end-1).^2));

            otherwise
             disp('curve not found, try another')
     end    
        % if len1 ~=0
        %     x = x/len*len1; dx = dx/len*len1; d2x = d2x/len*len1; len = len1;
        % end    
        nx = [dx(2,:);-dx(1,:)];
        nx = nx./sqrt(dx(1,:).^2+dx(2,:).^2);
        kappa = (dx(1,:).*d2x(2,:)-d2x(1,:).*dx(2,:))./(dx(1,:).^2+dx(2,:).^2).^1.5;
        len = 0; dom_area = 0;

        for ps = 1:pieces
            t = ts{ps}; t1 = ts1{ps}; dw = Dws{ps};    
            dw(isnan(dw)) = 0;
            dw((end+1)/2+1:end) = dw((end+1)/2-1:-1:1);
            t((end+1)/2+1:end) = -t((end+1)/2-1:-1:1)+2*pi;
            t1((end+1)/2+1:end) = -t1((end+1)/2-1:-1:1)+2*pi;
            ts{ps} = t(1:end-1); ts1{ps} = t1(1:end-1); Dws{ps} = dw(1:end-1);
            len = len + 2*pi/(length(ts{ps}))*sum(sqrt(dx(1,indices(ps,1):indices(ps,2)).^2+dx(2,indices(ps,1):indices(ps,2)).^2).*Dws{ps});
            dom_area = dom_area + 2*pi/length(ts{ps})*sum((dx(1,indices(ps,1):indices(ps,2)).*x(2,indices(ps,1):indices(ps,2))).*Dws{ps});
        end 
        % Dw = Dw(1:end-1); t1 = t1(1:end-1); t = t(1:end-1); 
        dom_area = abs(dom_area);
end