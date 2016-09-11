clear;close all;

% Set the number of grid points and build a cell-center grid
N=input(' Enter N, cell number - ');
L=10;
hb=1;
m=1;
p=2*pi;
sigma=2;
h=2*L/(N-1);
x=-L:h:L;
x=x';  % Turn x into a column vector.

%  Load the V coefficient array (make it a column vector)
V=zeros(size(x)); % (just 0 for now--we'll change it later)

% Load Dm with average values D(j-1/2) and Dp with D(j+1/2)
% Dm=zeros(N+2,1);Dp=zeros(N+2,1);  % Make the column vectors
% Dm(2:N+1)=.5*(D(2:N+1)+D(1:N));   % average j and j-1
% Dp(2:N+1)=.5*(D(2:N+1)+D(3:N+2)); % average j and j+1

% Set the initial distribution
Psi=1./(sqrt(sigma.*sqrt(pi))).*exp(1i.*p.*x./hb).*exp(-x.^2./(2.*sigma.^2));

% Find the maximum of T for setting plot limits
% Psimax=max(Psi);
% Psimin=min(Psi);
Psimin=-1;
Psimax=1;

% Choose the time step tau.
% The max tau for explicit stability is a reasonable choice
tau = input(' Enter the time step - ');

% Create the matrices A and B by loading them with zeros
A=zeros(N);
B=zeros(N);

% load A and B at interior points
for j=2:N-1
   A(j,j-1)= -1;
   A(j,j)  = -h^2*4*m/hb^2*1i*hb^2/tau+4*m*h^2/hb^2*V(j)/2+2;
   A(j,j+1)= -1;

   B(j,j-1)= 1;
   B(j,j)  = -h^2*4*m/hb^2*1i*hb^2/tau-4*m*h^2/hb^2*V(j)/2-2;
   B(j,j+1)= 1;
   
%    A(j,j-1)= hb^2/(4*m*h^2);
%    A(j,j)  = 1i*hb/tau-hb^2/(4*m)*2/h^2-V(j)/2;
%    A(j,j+1)= hb^2/(4*m*h^2);
% 
%    B(j,j-1)= -hb^2/(4*m*h^2);
%    B(j,j)  = 1i*hb/tau+hb^2/(4*m)*2/h^2+V(j)/2;
%    B(j,j+1)= -hb^2/(4*m*h^2);
end

% load the boundary conditions into A and B
%A(1,1)=0.5; A(1,2)=0.5; B(1,1)=0.;            % T(0)=0
%A(N+2,N+1)=0.5; A(N+2,N+2)=0.5; B(N+2,N+2)=0; % T(L)=0

% load the boundary conditions into A and B
A(1,1)=1; B(1,1)=0; % T(0)=0
A(N,N)=1; B(N,N)=0; % T(L)=0

% Set the number of time steps to take.
tfinal=input(' Enter the total run time - ');
nsteps=tfinal/tau;

% Choose how many iterations to skip between plots
skip = input(' Number of iterations to skip between plots - ');
times=zeros(nsteps,1);
expects=zeros(nsteps,1);
% This is the time advance loop.
for mtime=1:nsteps
    % define the time
    t=mtime*tau;
    % find the right-hand side for the solve at interior points
    r=B*Psi;

    % apply the boundary conditions
    r(1)=0;   % T(0)=0
    r(N)=0; % T(L)=0

    % do the linear solve to update T
    Psi=A\r;

    % Make a plot of T every once in a while.
    if(rem(mtime,skip) == 0)
        plot(x,real(Psi),'r-',x,abs(Psi).^2,'b-');
        % title(sprintf('Schrodinger Solution. Integral: %f',(trapz(abs(Psi).^2))*h))
        title(sprintf('Time: %f Integral: %f  Expectation: %f',t,(trapz(abs(Psi).^2))*h,(trapz(x.*abs(Psi).^2))*h))
        legend('real(Psi)','Psi^2');
        axis([-L L Psimin Psimax])
        pause(.1)
    end
    times(mtime)=t;
    expects(mtime)=(trapz(x.*abs(Psi).^2))*h;
end
figure;
plot(times,expects);
title('Expectation values over time')
axis([0 max(times) -L L])

