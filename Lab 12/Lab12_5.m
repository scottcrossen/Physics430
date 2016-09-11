% Staggered Leapfrog Script Template
clear;close all;

% Set the values for parameters
v0=1;

% build a cell-centered grid with N=200 cells
% on the interval x=0 to x=L, with L=1
N = 400;
xmin = 0;
xmax = 10;
L = xmax;
h = (xmax-xmin)/N; % spatial step
x = xmin-h/2:h:xmax+h/2;
j = 2:N+1; % spatial array to express indecies


% define the initial displacement and velocity vs. x
rho = 1+exp(-200*(x/L-1/2).^2);
% for problem 6 f
% rho = zeros(length(x),1);
% for i = 1:length(rho)
%     if x(i) <= L/2
%         rho(i) = 1.0;
%     end
% end
v = v0;

% Choose a time step
taulim=h/v0;
fprintf(' Courant time step limit %g \n',taulim)
tau=input(' Enter the time step - ')

% Get the initial value of yold from the initial y and vy
rhoold=rho;

% Apply the boundary conditions for yold(1) and yold(N+2)
rhoold(1)=2-rhoold(2);
rhoold(N+2)=2*rhoold(N+1)-rhoold(N); 

% plot the initial conditions and pause to look at them
subplot(2,1,1)
plot(x,rho)
xlabel('x');ylabel('y(x,0)');title('Initial Displacement')
subplot(2,1,2)
plot(x,v)
xlabel('x');ylabel('v_y(x,0)');title('Initial Velocity')
pause;

% Choose how long to run and when to plot
tfinal=input(' Enter tfinal - ')
skip=input(' Enter # of steps to skip between plots (faster) - ')
nsteps=tfinal/tau;

% here is the loop that steps the solution along

ymax = 1:nsteps;

figure  % open a new frame for the animation
for n=1:nsteps
   time=n*tau;  % compute the time

   % Use leapfrog and the boundary conditions to load ynew with y
   % at the next time step using y and yold, i.e., ynew(2:N+1)=...
   % Be sure to use colon commands so it will run fast.
   rhonew(j)=rho(j)-tau*v0/2/h*(rho(j+1)-rho(j-1))+v0^2*tau^2/(2*h^2)*(rho(j+1)-2*rho(j)+rho(j-1));
   
   % boundary condition
   rhonew(N+2)=2*rhonew(N+1)-rhonew(N);
   rhonew(1)=2-rhonew(2);
   
   %update yold and y
   rhoold=rho;rho=rhonew;

% make plots every skip time steps
   if mod(n,skip)==0
      plot(x,rho,'b-')
      xlabel('x');ylabel('y');
      title(['Staggered Leapfrog Wave: time=' num2str(time)])
      axis([min(x) max(x) 0 2.5]);
      pause(.1)
   end
   ymax(n)=max(abs(rhonew));
end