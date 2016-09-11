L=1.54;
N=500;
xmin=-5; xmax=5;
h=(xmax-xmin)/(N);
x=xmin-h/2:h:xmax+h/2;
A=zeros(N+2,N+2);
A(1,1)=1/2;
A(1,2)=1/2;
A(N+2,N+1)=1/2;
A(N+2,N+2)=1/2;
B=eye(N+2,N+2);
B(1,1)=0;
B(N+2,N+2)=0;
for iterator1=2:1:N+1
    A(iterator1,iterator1-1)=-1/2*1/(h^2);
    A(iterator1,iterator1)=-1/2*-2/(h.^2)+x(iterator1)^4;
    A(iterator1,iterator1+1)=-1/2*1/(h^2);
end;

[V,D]=eig(A,B);  % find the eigenvectors and eigenvalues

w2raw=diag(D);  % convert lambda to omega^2

[w2,k]=sort(w2raw);  % sort omega^2 into ascending along with a
                     % sort key k(n) that remembers where each
                     % omega^2 came from so we can plot the proper
                     % eigenvector in V

for n=1:N  % run through the sorted list and plot each eigenvector
   % load the plot title into t
   t=sprintf(' epsilon = %g',w2(n));
   gn=V(:,k(n)); % extract the eigenvector
   plot(x,abs(gn).^2,'b-');  % plot the eigenvector that goes with omega^2
   title(t);xlabel('x');ylabel('g(n,x)');  % label the graph
   pause
end


