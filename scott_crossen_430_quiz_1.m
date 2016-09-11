% Name: Scott Leland Crossen
% Take Home Problem #1
% Physics 430 Section 2

% Progam Output:
% Numeric k for  n=1 is 2.404756. 	Abs(Error) in k_1 is 0.000029
% Numeric k for  n=2 is 5.519340. 	Abs(Error) in k_2 is 0.000134
% Numeric k for  n=3 is 8.650940. 	Abs(Error) in k_3 is 0.000322
% Numeric k for  n=4 is 11.784527. 	Abs(Error) in k_4 is 0.000594

clear; close all;
N=100; % Amount of points
% This next line was found using mathematica and are the first four zeros of J0[x].
kactual=[2.404825557695773; 5.5200781102863115; 8.653727912911013; 11.791534439014281];
xmin=0; xmax=1; % Plot between zero and one.
h=(xmax-xmin)/(N-1);
x=xmin:h:xmax; % Set domain array
A=zeros(N,N);
% Set derivative boundary condition:
A(1,1)=3/(2*h);
A(1,2)=-2/h;
A(1,3)=1/(2*h);
% Set value boundary condition:
A(N,N)=1;
% Create 'B' Matrix:
B=eye(N,N);
B(1,1)=0;
B(N,N)=0;
% Map the differntial equation onto 'A':
for iterator1=2:1:N-1
    A(iterator1,iterator1-1)=1/(h^2)-1/(2*h)/x(iterator1);
    A(iterator1,iterator1)=-2/(h.^2);
    A(iterator1,iterator1+1)=1/(h^2)+1/(2*h)/x(iterator1);
end;
[V,D]=eig(A,B);
k2raw=-diag(D);
[k2,otherk]=sort(k2raw); % Sort Eigenvalues^2.
for n=3:1:6 % Plot the first four eigenvalues.
figure;
gn=V(:,otherk(n));
t=sprintf('k = %g ',sqrt(abs(k2(n)))); % Display the eigenvalue as a title.
plot(x,gn./gn(1),'b.',x,besselj(0,sqrt(abs(k2(n))).*x),'r-'); % Plot with the actual result.
% Display the error of the calculated eigenvalue with the one found in mathematica.
fprintf('Numeric k for  n=%d is %f. \tAbs(Error) in k_%d is %f\n',n-2,sqrt(abs(k2(n))),n-2,abs(sqrt(abs(k2(n)))-kactual(n-2))/kactual(n-2));
title(t);xlabel('x');ylabel('y(x)');
legend('Approximate','Actual');
end;