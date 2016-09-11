clear; close all;
N=500;
a=0;b=pi/2;
h=(b-a)/(N);
x=a-h/2:h:b+h/2;
f=sin(x);
plot(x,f);
% a function-value BC can be imposed by setting the two end point's average
% value to the specified condition
% a derivative-value BC can be imposed by setting the slope between the two
% lines to the specified value. (derivative =0 -> points equal)
