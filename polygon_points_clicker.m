figure
plot(0,0)
xlim([0 5]);ylim([0 5]); grid on
[x,y] = getpts;%input("whats the point yo");
curve_params = [x y]; 
writematrix(curve_params,"poly_points.csv")