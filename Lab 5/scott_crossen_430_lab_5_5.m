% Staggered Leapfrog Script Template
clear;close all;

% Set the values for parameters
L=1.2;
T=127;
xmin=0;
xmax=L;
gam=50;
mu=.003;
omega=1080;
c=sqrt(T/mu); % wave speed

% build a cell-centered grid with N=200 cells
% on the interval x=0 to x=L, with L=1
N=200;
J=2:N+1;
h=(xmax-xmin)/(N);
x=xmin-h/2:h:xmax+h/2;
f=zeros(size(x));
for tempx=1:length(x)
    if (.8 <= x(tempx) && x(tempx) <=1)
    f(tempx)=.73;
    end;
end;

% define the initial displacement and velocity vs. x
y = zeros(1,length(x));
vy = y;

% Choose a time step (Suggest that it is no larger than taulim)
taulim=h/c;
fprintf(' Courant time step limit %g \n',taulim);
tau=input(' Enter the time step - ');

% Get the initial value of yold from the initial y and vy
%yold(J)=(2+gam*tau)/(2*gam*tau)*(-2*tau*vy(J)+(1/(2+gam*tau)*(4*y(J)+2*c^2*tau^2/(h^2)*(y(J+1)-2*y(J)+y(J-1)))));
yold(J)=0*J;

% Apply the boundary conditions for yold(1) and yold(N+2)
yold(1)=-yold(2);
yold(N+2)=-yold(N+1);

% plot the initial conditions and pause to look at them
subplot(2,1,1)
plot(x,y)
xlabel('x');ylabel('y(x,0)');title('Initial Displacement')
subplot(2,1,2)
plot(x,vy)
xlabel('x');ylabel('v_y(x,0)');title('Initial Velocity')
pause;

% Choose how long to run and when to plot
tfinal=input(' Enter tfinal - ');
skip=input(' Enter # of steps to skip between plots (faster) - ');
nsteps=tfinal/tau;

% here is the loop that steps the solution along
times=(1:nsteps)*tau;
%times=times*tau;
amps=zeros(size(times));
figure  % open a new frame for the animation
for n=1:nsteps
   time=n*tau;  % compute the time

   % Use leapfrog and the boundary conditions to load ynew with y
   % at the next time step using y and yold, i.e., ynew(2:N+1)=...
   % Be sure to use colon commands so it will run fast.
   % ynew(J)=1/(2+gam*tau)*(4*y(J)-2*yold(J)+gam*tau*yold(J)+2*c^2*tau^2/(h^2)*(y(J+1)-2*y(J)+y(J-1)));
   ynew(J)=1/(2+gam*tau)*(4*y(J)-2*yold(J)+gam*tau*yold(J)+2*c^2*tau^2/(h^2)*(y(J+1)-2*y(J)+y(J-1))+2*f(J)/mu*tau^2*cos(time*omega));
%ynew(J)=(mu*(2*c^2*y(J+1)-2*y(J)+y(J-1)*tau^2+h^2*(4*y(J)-2*yold(J)+yold(J)*gam*tau))+2*f(J)*h^2*tau^2*cos(time*omega))/(h^2*mu*(2+gam*tau));
ynew(1)=-ynew(2);
ynew(N+2)=-ynew(N+1);

   %update yold and y
   yold=y;y=ynew;

% make plots every skip time steps
   if mod(n,skip)==0
      plot(x,y,'b-')
      xlabel('x');ylabel('y');
      title(['Staggered Leapfrog Wave: time=' num2str(time)])
      axis([min(x) max(x) -5*10^-4 5*10^-4]);
      pause(.1)
   end
   amps(n)=max(abs(ynew));
end
%plot(times,amps,'b-');
%title('Max Amplitude vs. Time');
