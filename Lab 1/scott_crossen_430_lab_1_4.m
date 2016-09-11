clear; close all;
N=100;
a=0;b=5;
h=(b-a)/(N-1);
x=a:h:b;
f=besselj(0,x);
fp(2:N-1)=(f(3:N)-f(1:N-2))/(2*h);
fpp(2:N-1)=(f(3:N)-2*f(2:N-1)+f(1:N-2))/h^2;

% Linear Extrapolation
fp(1)=2*fp(2)-fp(3);
fp(N)=2*fp(N-1)-fp(N-2);
fpp(1)=2*fpp(2)-fpp(3);
fpp(N)=2*fpp(N-1)-fpp(N-2);
plot(x,f,x,fp,x,fpp,x,-besselj(1,x),x,.5*(-besselj(0,x)+besselj(2,x)));
title('Linear Extrapolation');
legend('Function','Approx Deriv.','Approx 2 Deriv.','Exact Deriv.','Exact 2 Deriv.')
fprintf('For Linear Extrapolation\n');
fprintf('at f(0) df/dx=%f and d^2f/dx^2=%f\n',fp(1),fpp(1));
fprintf('at f(5) df/dx=%f and d^2f/dx^2=%f\n',fp(N),fpp(N));

% Quadratic Extrapolation
figure;
fp(1)=3*fp(2)-3*fp(3)+fp(4);
fp(N)=3*fp(N-1)-3*fp(N-2)+fp(N-3);
fpp(1)=3*fpp(2)-3*fpp(3)+fpp(4);
fpp(N)=3*fpp(N-1)-3*fpp(N-2)+fpp(N-3);
plot(x,f,x,fp,x,fpp,x,-besselj(1,x),x,.5*(-besselj(0,x)+besselj(2,x)));
title('Quadratic Extrapolation');
legend('Function','Approx Deriv.','Approx 2 Deriv.','Exact Deriv.','Exact 2 Deriv.')
fprintf('For Quadratic Extrapolation\n');
fprintf('at f(0) df/dx=%f and d^2f/dx^2=%f\n',fp(1),fpp(1));
fprintf('at f(5) df/dx=%f and d^2f/dx^2=%f\n',fp(N),fpp(N));

% Actual Values at Boundary Conditions
fprintf('For the Actual Values\n');
fprintf('at f(0) df/dx=0 and d^2f/dx^2=.327579\n');
fprintf('at f(5) df/dx=-.5 and d^2f/dx^2=.112081\n');
