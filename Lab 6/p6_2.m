% Staggered Leapfrog Script Template
clear;close all;

% Set the values for parameters
sigma=2;
mu=.3;
c=sqrt(sigma/mu);

% build a cell-centered grid with N=200 cells
N=100;
j = 2:N-1; % spatial array to express internal indecies
xmin=-5; xmax=5;
hx=(xmax-xmin)/(N-1);
x=xmin:hx:xmax;

ymin=-5; ymax=5;
hy=(ymax-ymin)/(N-1);
y=ymin:hy:ymax;

h=hx;
% make 2d spatial grid
[X,Y]=ndgrid(x,y);


% define the initial displacement and velocity vs. z
% for parts a-c
% z=exp(-5*(X.^2+Y.^2));
% vz=zeros(length(x),length(y));
% for part d
vz=exp(-5*(X.^2+Y.^2));
z=zeros(length(x),length(y));
% pin down edges to zero
z(:,1)=0;
z(:,N)=0;
z(1,:)=0;
z(N,:)=0;

% Choose a time step (Suggest that it is no larger than taulim)
% here there is a factor of sqrt(2) in front of h/c
taulim=.69*h/c;
fprintf(' Courant time step limit %g \n',taulim);
tau=input(' Enter the time step - ');

% Get the initial value of zold from the initial z and vz
zold(j,j) = -vz(j,j)*tau + z(j,j) + tau^2*sigma*((z(j+1,j)-2*z(j,j)+z(j-1,j)) + (z(j,j+1)-2*z(j,j)+z(j,j-1)))/(2*mu*h^2);

% Apply the boundary conditions for yold(1) and yold(N+2)
zold(:,1)=0;
zold(:,N)=0;
zold(1,:)=0;
zold(N,:)=0;

figure;
% plot the initial conditions and pause to look at them
subplot(2,1,1);
surf(x,y,z);
xlabel('x');ylabel('y(x,0)');title('Initial Displacement');
subplot(2,1,2);
surf(x,y,vz);
xlabel('x');ylabel('v_y(x,0)');title('Initial Velocity');
pause;

% Choose how long to run and when to plot
tfinal=input(' Enter tfinal - ');
skip=input(' Enter # of steps to skip between plots (faster) - ');
nsteps=tfinal/tau;

% here is the loop that steps the solution along

% define an array for z(0,0,t)
zZeros=zeros(1,floor(nsteps));
figure  % open a new frame for the animation
for n=1:nsteps
   time=n*tau;  % compute the time

   % Use leapfrog and the boundary conditions to load ynew with y
   % at the next time step using y and yold, i.e., ynew(2:N+1)=...
   % Be sure to use colon commands so it will run fast.
   znew(j,j)=2*z(j,j) - zold(j,j) + tau^2*sigma*((z(j+1,j)-2*z(j,j)+z(j-1,j)) + (z(j,j+1)-2*z(j,j)+z(j,j-1)))/(mu*h^2);
   
   % boundary conditions
   znew(:,1)=0;
   znew(:,N)=0;
   znew(1,:)=0;
   znew(N,:)=0;
   %update yold and y
   zold=z;z=znew;
   
   zZeros(1,n)=znew(floor(N/2),floor(N/2));

% make plots every skip time steps
   if mod(n,skip)==0
      surf(x,y,z)
      xlabel('x');ylabel('y');
      title(['Staggered Leapfrog Wave: time=' num2str(time)])
      caxis([-.25 .25])
      pause(.1)
   end
end

% plot z(0,0,t)
figure
times=tau:tau:tfinal;
plot(times,zZeros);
title('Z(0,0,t)')