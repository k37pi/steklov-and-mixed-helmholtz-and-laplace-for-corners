function [x, dx, d2x] = slantline(t,spx,epx,spy,epy) 
    mx = (epx-spx)/(2*pi); my = (epy-spy)/(2*pi); 
    x = [spx+mx*t;spy+my*t];
    mp = x(1,(end+1)/2);
    x(1,((end+1)/2+1):end) = 2*mp-x(1,((end+1)/2-1):-1:1);
    mp = x(2,(end+1)/2);
    x(2,((end+1)/2+1):end) = 2*mp-x(2,((end+1)/2-1):-1:1);
    dx = [mx+0*t;my+0*t];
    d2x = [0*t;0*t];
    x = x(:,1:end-1); dx = dx(:,1:end-1); d2x = d2x(:,1:end-1);
end


           