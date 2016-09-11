Nx=30;
a=0; b=2;
dx=(b-a)/(Nx-1);
x=a:dx:b;

Ny=50;
c=-1; d=3;
dy=(d-c)/(Ny-1);
y=c:dy:d;

[X,Y]=ndgrid(x,y);

f=exp(-(X.^2+Y.^2)).*cos(5*sqrt(X.^2.+Y.^2));

% plots by index
surf(f)
title('surf(f)')
figure
% proide grid for it to plot on
surf(Y,X,f)
title('surf(Y,X,f)')
figure
% it constructs a grid to plot on
surf(y,x,f)
title('surf(y,x,f)')