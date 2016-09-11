clear; close all;
N=100;
a=0;b=pi;
h=(b-a)/(N-1);
x=a:h:b;
f=sin(x).*sinh(x);
plot(x,f);
whos('x*')