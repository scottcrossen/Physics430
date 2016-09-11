clear; close all;
N=5000;
a=0;b=2;
h=(b-a)/(N);
x=a+h/2:h:b-h/2;
f=cos(x);
plot(x,f);
fprintf('the approximate area is about %f\n',sum(f)*h);
fprintf('the actual area is about %f\n',0.909297);
