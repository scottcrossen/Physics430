% Staggered Leapfrog Script Template
clear;close all;

% Set the values for parameters
D=4;
L=3;
C=1/2;
xmin=0;
xmax=L;
gam=0.2;

% build a cell-centered grid with N=200 cells
% on the interval x=0 to x=L, with L=1
N=80;
J=2:N+1;
h=(xmax-xmin)/(N);
x=xmin-h/2:h:xmax+h/2;

% define the initial displacement and velocity vs. x
%T = sin(pi*x/L);
T = exp(-40*(x/L-1/2).^2);
vT = 0;

% Choose a time step (Suggest that it is no larger than taulim)
taulim=C*h^2/D;
fprintf(' Courant time step limit %g \n',taulim);
tau=input(' Enter the time step - ');

% Apply the boundary conditions for T(1) and T(N+2)
% The first set of BCs is for value BCs. the second for derivative BCs.
T(1) = -T(2);
T(N+2) = -T(N+1);
% T(1) = T(2);
% T(N+2) = T(N+1);

% plot the initial conditions and pause to look at them
subplot(2,1,1)
plot(x,T)
xlabel('x');ylabel('T(x,0)');title('Initial Displacement')
subplot(2,1,2)
plot(x,vT)
xlabel('x');ylabel('T_y(x,0)');title('Initial Velocity')
pause;

% Choose how long to run and when to plot
tfinal=input(' Enter tfinal - ');
skip=input(' Enter # of steps to skip between plots (faster) - ');
nsteps=tfinal/tau;

% here is the loop that steps the solution along
times=(1:nsteps)*tau;
%times=times*tau;
errors=zeros(size(times));
figure  % open a new frame for the animation
for n=1:nsteps
   time=n*tau;  % compute the time

   % Use leapfrog and the boundary conditions to load ynew with y
   % at the next time step using y and yold, i.e., ynew(2:N+1)=...
   % Be sure to use colon commands so it will run fast.
    Tnew(J)=T(J)+D*tau/h^2*(T(J+1)-2*T(J)+T(J-1));
% The first set of BCs is for value BCs. the second for derivative BCs.
    Tnew(1) = -Tnew(2);
    Tnew(N+2) = -Tnew(N+1);
%     Tnew(1) = Tnew(2);
%     Tnew(N+2) = Tnew(N+1);

   %update yold and y
   T=Tnew;

% make plots every skip time steps
   if mod(n,skip)==0
      plot(x,T,'r.')
      xlabel('x');ylabel('y');
      title(['Staggered Leapfrog Wave: time=' num2str(time)])
      axis([min(x) max(x) -2 2]);
      pause(.1)
   end
end
