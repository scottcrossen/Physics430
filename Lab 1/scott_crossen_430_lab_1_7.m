clear; close all;
N=100;
a=0;b=5;
h=(b-a)/(N-1);
x=a:h:b;
f=cos(x)+.001*rand(1,length(x));
fp(2:N-1)=(f(3:N)-f(1:N-2))/(2*h);
fpp(2:N-1)=(f(3:N)-2*f(2:N-1)+f(1:N-2))/h^2;
fp(1)=3*fp(2)-3*fp(3)+fp(4);
fp(N)=3*fp(N-1)-3*fp(N-2)+fp(N-3);
fpp(1)=3*fpp(2)-3*fpp(3)+fpp(4);
fpp(N)=3*fpp(N-1)-3*fpp(N-2)+fpp(N-3);

gp=-sin(x);
gpp=-cos(x);
plot(x,f,x,fp);%,x,fpp)
legend('Function','First Derivative','Second Derivative');
title('Function and Approximate Derivatives');


figure;
plot(x,fp-gp);%,x,fpp-gpp);
legend('First Derivative','Second Derivative');
title('Difference Between Approximate and Actual Derivatives');