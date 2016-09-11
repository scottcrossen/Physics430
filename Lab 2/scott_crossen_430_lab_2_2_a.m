clear; close all;

N=30;
xmin=0; xmax=5;
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax;
x=x';
y1=-4./besselj(1,5).*besselj(1,x)+x;

A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;
b=zeros(N,1);
b(1,1)=0;
b(N,1)=1;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2)-1/(2*h.*x(iterator1));
    A(iterator1,iterator1)=-2/(h.^2)+1-1./(x(iterator1).^2);
    A(iterator1,iterator1+1)=1/(h^2)+1/(2*h.*x(iterator1));
    b(iterator1,1)=x(iterator1);
end;
y2=A\b;



figure;
plot(x,y1,x,y2,'r.');
legend('Exact','Approximation');
title('Equation 2.6');

figure;
plot(x,y1-y2);
title('Equation 2.6 Error');
legend('Exact-Approx');