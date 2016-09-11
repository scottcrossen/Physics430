mu=.003;
T=127;
L=1.2;
N=500;
xmin=0; xmax=1.2;
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
for n=1:N
   t=sprintf(' w^2 = %g w = %g ',w2(n),sqrt(abs(w2(n))) );
   gn=V(:,k(n));
   plot(x,gn,'b-');
   title(t);xlabel('x');ylabel('g(n,x)');
   pause
end

