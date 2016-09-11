clear; close all;

N=30;
xmin=0; xmax=3;
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax;
x=x';
y1=.95*(x-1.5).^2-2.13;
y2=y1;

A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;
b=zeros(N,1);
b(1,1)=0;
b(N,1)=0;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2);
    A(iterator1,iterator1)=-2/(h.^2);
    A(iterator1,iterator1+1)=1/(h^2);
    b(iterator1,1)=1-sin(y2(iterator1));
end;

err=1;
while err>.0001
    y2=A\b;
    err1=A*y2-(1-sin(y2));
    err=max(abs(err1(2:N-1)));
    b(2:N-1)=1-sin(y2(2:N-1));
end;

figure;
plot(x,y1,x,y2,'r.');
legend('Exact','Approximation');
title('Equation 2.13');

figure;
plot(x,y1-y2);
title('Equation 2.13 Error');
legend('Exact-Approx');