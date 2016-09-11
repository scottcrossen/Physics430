clear; close all;

N=30;
xmin=0; xmax=2;
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax;
x=x';
y1=x./9-sin(3.*x)/(27.*cos(6));

A=zeros(N,N);
A(1,1)=1;
A(N,N-1)=-1;
A(N,N)=1;
b=zeros(N,1);
b(1,1)=0;
b(N,1)=0;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2);
    A(iterator1,iterator1)=9-2/(h.^2);
    A(iterator1,iterator1+1)=1/(h^2);
    b(iterator1,1)=x(iterator1);
end;
y2=A\b;


A=zeros(N,N);
A(1,1)=1;
A(N,N-2)=1/(2*h);
A(N,N-1)=-2/h;
A(N,N)=3/(2*h);
b=zeros(N,1);
b(1,1)=0;
b(N,1)=0;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2);
    A(iterator1,iterator1)=9-2/(h.^2);
    A(iterator1,iterator1+1)=1/(h^2);
    b(iterator1,1)=x(iterator1);
end;
y3=A\b;



figure;
plot(x,y1,x,y3,'r.');
legend('Exact','Approximation');
title('Equation 2.10');

figure;
plot(x,y1-y2,x,y1-y3);
title('Equation 2.10 Error');
legend('Exact-Linear End', 'Exact-Quad End');