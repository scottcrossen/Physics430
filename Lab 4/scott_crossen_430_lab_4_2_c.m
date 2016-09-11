L=1.54;
N=500;
xmin=0; xmax=1.54;
h=(xmax-xmin)/(N);
x=xmin-h/2:h:xmax+h/2;
A=zeros(N+2,N+2);
A(1,1)=-1/h;
A(1,2)=1/h;
A(N+2,N+1)=1/2;
A(N+2,N+2)=1/2;
B=eye(N+2,N+2);
B(1,1)=1/2;
B(1,2)=1/2;
B(N+2,N+2)=0;
for iterator1=2:1:N+1
    A(iterator1,iterator1-1)=x(iterator1)/(h^2)-1/(2*h);
    A(iterator1,iterator1)=-2*x(iterator1)/(h.^2);
    A(iterator1,iterator1+1)=x(iterator1)/(h^2)+1/(2*h);
end;

[V,D]=eig(A,B);  % find the eigenvectors and eigenvalues

g=9.8;

w2raw=-(g)*diag(D);  % convert lambda to omega^2

[w2,k]=sort(w2raw);  % sort omega^2 into ascending along with a
                     % sort key k(n) that remembers where each
                     % omega^2 came from so we can plot the proper
                     % eigenvector in V

% Change start to 2 from 1 for start.
for n=2:N  % run through the sorted list and plot each eigenvector
   % load the plot title into t
   t=sprintf(' w^2 = %g w = %g ',w2(n),sqrt(abs(w2(n))) );
   gn=V(:,k(n)); % extract the eigenvector
   plot(x,gn,'b-');  % plot the eigenvector that goes with omega^2
   title(t);xlabel('x');ylabel('g(n,x)');  % label the graph
   pause
end


