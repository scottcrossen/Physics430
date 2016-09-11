mu=.003;
T=127;
L=1.2;
N=500;
xmin=0; xmax=L;
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax;
A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;
B=eye(N,N);
B(1,1)=0;
B(N,N)=0;
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2);
    A(iterator1,iterator1)=-2*1/(h.^2);
    A(iterator1,iterator1+1)=1/(h^2);
end;
[V,D]=eig(A,B);
w2raw=-(T/mu)*diag(D);
[w2,k]=sort(w2raw);

n=3;
figure;
t=sprintf(' w^2 = %g w = %g ',w2(n),sqrt(abs(w2(n))) );
gn=V(:,k(n));
plot(x,gn,'b-',x,max(gn)*sin((n-2)*pi/L*sqrt(T/mu)*sqrt(mu/T)*x),'r-');
title(t);xlabel('x');ylabel('g(n,x)');
legend('Approximate','Actual');
fprintf('Actual  Omega for  n=%d is %f\n',n-2,(n-2)*pi/L*sqrt(T/mu));
fprintf('Numeric Omega for  n=%d is %f\n',n-2,sqrt(w2(n)));

n=4;
figure;
t=sprintf(' w^2 = %g w = %g ',w2(n),sqrt(abs(w2(n))) );
gn=V(:,k(n));
plot(x,gn,'b-',x,max(gn)*sin((n-2)*pi/L*sqrt(T/mu)*sqrt(mu/T)*x),'r-');
legend('Approximate','Actual');
title(t);xlabel('x');ylabel('g(n,x)');
fprintf('Actual  Omega for  n=%d is %f\n',n-2,(n-2)*pi/L*sqrt(T/mu));
fprintf('Numeric Omega for  n=%d is %f\n',n-2,sqrt(w2(n)));

n=5;
figure;
t=sprintf(' w^2 = %g w = %g ',w2(n),sqrt(abs(w2(n))) );
gn=V(:,k(n));
plot(x,gn,'b-',x,max(gn)*sin((n-2)*pi/L*sqrt(T/mu)*sqrt(mu/T)*x),'r-');
legend('Approximate','Actual');
title(t);xlabel('x');ylabel('g(n,x)');
fprintf('Actual  Omega for  n=%d is %f\n',n-2,(n-2)*pi/L*sqrt(T/mu));
fprintf('Numeric Omega for  n=%d is %f\n',n-2,sqrt(w2(n)));


n=22;
figure;
t=sprintf(' w^2 = %g w = %g ',w2(n),sqrt(abs(w2(n))) );
gn=V(:,k(n));
plot(x,gn,'b-',x,max(gn)*sin((n-2)*pi/L*sqrt(T/mu)*sqrt(mu/T)*x),'r-');
legend('Approximate','Actual');
title(t);xlabel('x');ylabel('g(n,x)');
fprintf('Actual  Omega for  n=%d is %f\n',n-2,(n-2)*pi/L*sqrt(T/mu));
fprintf('Numeric Omega for  n=%d is %f\n',n-2,sqrt(w2(n)));
