clear; close all;

N=500;
xmin=0; xmax=1.2;
mu=.003;
T=127;
L=1.2;
w=400;
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax;

A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;
b=zeros(N,1);
b(1,1)=0;
b(N,1)=0;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=T/(h^2);
    A(iterator1,iterator1)=-2*T/(h.^2)+mu*w^2;
    A(iterator1,iterator1+1)=T/(h^2);
    if (.8 <= x(iterator1)) && (x(iterator1) <=1)
        b(iterator1,1)=-.73;
    end;
end;
y2=A\b;
figure;
plot(x,y2,'r.');
title('g(x)');

figure;
omegas=400:(1200-400)/100:1200;
m=zeros(size(omegas));
for iterator2=1:101
    for iterator1=2:1:N-1
        A(iterator1,iterator1-1)=T/(h^2);
        A(iterator1,iterator1)=-2*T/(h.^2)+mu*omegas(iterator2)^2;
        A(iterator1,iterator1+1)=T/(h^2);
    end;
    y2=A\b;
    plot(x,y2,'r-');
    title('g(x)');
    axis([0 1.4 -4.5*10^-4 4.5*10^-4])
    pause(.05);
    m(iterator2)=max(abs(y2));
end;
figure;
plot(omegas,m,'r-');
title('Max. Amplitude');
fprintf('Relative amplitudes at Omega = %f and Omega= %f\n',max(m),max(m(50:end)))