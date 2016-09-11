clear; close all;

Nx=30;
a=0;
b=2;
dx=(b-a)/(Nx-1);
x=a:dx:b;
Ny=50;
c=-1;
d=3;
dy=(d-c)/(Ny-1);
y=c:dy:d;
[X,Y]=meshgrid(x,y);

f=exp(-(X.^2+Y.^2)).*cos(5.*sqrt(X.^2+Y.^2));
surf(f);
xlabel('x');
ylabel('y');
figure;
surf(X,Y,f);
xlabel('x');
ylabel('y');
figure;
surf(x,y,f);
xlabel('x');
ylabel('y');