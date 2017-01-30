%---------------------------------------------------%
% AFIT Qualifying Exam Question #4, 01 Dec 2014     %
% Timothy Coon                                      %
% Advisor: Dr Cobb                                  %
%---------------------------------------------------%
% Q4_GradientPlot.m
% plot the "terrain" as defined by the (x,y)-dependent energy contribution
% in the energy cost functional.

syms x y
f = exp(-((x-2)^2+(y-2)^2));

gradf = jacobian(f,[x,y]);

[xx, yy] = meshgrid(x0:.1:xf,ymin:.1:ymax);
ffun = @(x,y) eval(vectorize(f));
fxfun = @(x,y) eval(vectorize(gradf(1)));
fyfun = @(x,y) eval(vectorize(gradf(2)));
contour(xx, yy, ffun(xx,yy), 30)
hold on
[xx, yy] = meshgrid(x0:.25:xf,ymin:.25:ymax);
quiver(xx, yy, fxfun(xx,yy), fyfun(xx,yy), 0.6)
axis equal tight, hold off