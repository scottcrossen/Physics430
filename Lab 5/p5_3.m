% Staggered Leapfrog Script Template
clear;close all;

% Set the values for parameters
c=2; % wave speed

% build a cell-centered grid with N=200 cells
% on the interval x=0 to x=L, with L=1
N = 200;
xmin = 0;
xmax = 1;
L = xmax;
h = (xmax-xmin)/N; % spatial step
x = xmin-h/2:h:xmax+h/2;
j = 2:N+1; % spatial array to express indecies



% define the initial displacement and velocity vs. x
%for parts a/b
%y = exp(-(x-L/2).^2*160/L^2)-exp(-(0-L/2).^2*160/L^2);
%vy = zeros(1,length(x));
%for part c
vy = exp(-(x-L/2).^2*160/L^2)-exp(-(0-L/2).^2*160/L^2);
y = zeros(1,length(x));

% Choose a time step (Suggest that it is no larger than taulim)
taulim=h/c;
fprintf(' Courant time step limit %g \n',taulim)
tau=input(' Enter the time step - ')

% Get the initial value of yold from the initial y and vy
yold(j) = y(j)-vy(j)*tau + c^2*tau^2*(y(j+1)-2*y(j)+y(j-1))/(2*h^2);

% Apply the boundary conditions for yold(1) and yold(N+2)
%for parts a/c
yold(1)=-yold(2);
yold(N+2)=-yold(N+1);
% for part b
%yold(1)=yold(2);
%yold(N+2)=yold(N+1);

% plot the initial conditions and pause to look at them
subplot(2,1,1)
plot(x,y)
xlabel('x');ylabel('y(x,0)');title('Initial Displacement')
subplot(2,1,2)
plot(x,vy)
xlabel('x');ylabel('v_y(x,0)');title('Initial Velocity')
pause;

% Choose how long to run and when to plot
tfinal=input(' Enter tfinal - ')
skip=input(' Enter # of steps to skip between plots (faster) - ')
nsteps=tfinal/tau;

% here is the loop that steps the solution along

figure  % open a new frame for the animation
for n=1:nsteps
   time=n*tau;  % compute the time

   % Use leapfrog and the boundary conditions to load ynew with y
   % at the next time step using y and yold, i.e., ynew(2:N+1)=...
   % Be sure to use colon commands so it will run fast.
   ynew(j)=2*y(j)-yold(j)+c^2*tau^2*(y(j+1)-2*y(j)+y(j-1))/h^2;
   % for parts a/c
   ynew(1)=-ynew(2);
   ynew(N+2)=-ynew(N+1);
   % for part b
   %ynew(1)=ynew(2);
   %ynew(N+2)=ynew(N+1);
   %update yold and y
   yold=y;y=ynew;

% make plots every skip time steps
   if mod(n,skip)==0
      plot(x,y,'b-')
      xlabel('x');ylabel('y');
      title(['Staggered Leapfrog Wave: time=' num2str(time)])
      % for parts a/b
      %axis([min(x) max(x) -2 2]);
      % for part c
      axis([min(x) max(x) -.1 .1]);
      pause(.1)
   end
end
