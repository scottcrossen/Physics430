clear; close all;

N=5000;
xmin=0; xmax=5;
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax;
x=x';

A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;
b=zeros(N,1);
b(1,1)=0;
b(N,1)=3;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2)-sin(x(iterator1))/(2*h);
    A(iterator1,iterator1)=-2/(h.^2)+exp(x(iterator1));
    A(iterator1,iterator1+1)=1/(h^2)+sin(x(iterator1))/(2*h);
    b(iterator1,1)=x(iterator1).^2;
end;
y2=A\b;



figure;
plot(x,y2,'r-');
legend('Approximation');
title('Equation 2.7');

fprintf('value at 4.5 is %f\n',y2(ceil(4.5./h)));