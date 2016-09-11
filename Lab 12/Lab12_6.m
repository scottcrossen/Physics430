clear;close all;

v0 = 1;

% Set the number of grid points and build a cell-center grid
N=input(' Enter N, cell number - ');
L=10;
hb=1;
m=1;
p=2*pi;
sigma=2;
h=L/(N);
x=-h/2:h:L+h/2;
x=x';  % Turn x into a column vector.

% Set the initial distribution
rho= 1+exp(-200.*(x/L-1/2).^2); % for parts a-d
% for part e
rho = zeros(length(x),1);
for i = 1:length(rho)
    if x(i) <= L/2
        rho(i) = 1.0;
    end
end

% v
v=v0*ones(length(rho),1); % for parts a-c, e-f
% v=1.2-x/L; % for part d

% Find the maximum of Rho for setting plot limits
rhomin=0;
rhomax=2.5;

% Choose the time step tau.
% The max tau for explicit stability is a reasonable choice
tau = input(' Enter the time step - ');

% Create the matrices A and B by loading them with zeros
A=zeros(N+2);
B=zeros(N+2);

% load A and B at interior points
for j=2:N+1
   A(j,j-1)= v(j-1);
   A(j,j)  = -4*h/tau;
   A(j,j+1)= -v(j+1);

   B(j,j-1)= -v(j-1);
   B(j,j)  = -4*h/tau;
   B(j,j+1)= v(j+1);
end

% load the boundary conditions into A and B
A(1,1)=0.5; A(1,2)=0.5; B(1,1)=0; % rho(0)=0
A(N+2,N+2)=1; A(N+2,N+1)=-2; A(N+2,N)=1; 
B(N+2,N+2)=0; % T(L)=0

% Set the number of time steps to take.
tfinal=input(' Enter the total run time - ');
nsteps=tfinal/tau;

% Choose how many iterations to skip between plots
skip = input(' Number of iterations to skip between plots - ');
times=zeros(nsteps,1);
% This is the time advance loop.
for mtime=1:nsteps
    % define the time
    t=mtime*tau;
    % find the right-hand side for the solve at interior points
    r=B*rho;

    % apply the boundary conditions
    r(1)=1;   % rho(0)=0
    r(N+2)=0; % rho(L)=0

    % do the linear solve to update rho
    rho=A\r;

    % Make a plot of rho every once in a while.
    if(rem(mtime,skip) == 0)
        plot(x,rho,'r-');
        axis([0 L rhomin rhomax])
        pause(.1)
    end
end

