clear; close all;

min=0;
max=10;
N=100;
tau=(max-min)/N;
gam=21;
%2/tau;
y=zeros(1,N);
time=zeros(1,N);
y(1)=1;
for n=2:N    
y(n)=y(n-1)*(1-gam*tau/2)/(1+gam*tau/2);
time(n)=n*tau;
end;
plot(time,y,time,exp(-gam*time));