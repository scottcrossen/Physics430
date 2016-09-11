% START OF PROBLEM 1
clear;close all;

% Set the values for parameters
sigma=2;
mu=.3;
c=sqrt(sigma/mu);
L=5;

% build a cell-centered grid with N=200 cells
N=100;
j = 2:N-1; % spatial array to express internal indecies
xmin=-L; xmax=L;
hx=(xmax-xmin)/(N-1);
x=xmin:hx:xmax;

ymin=-L; ymax=L;
hy=(ymax-ymin)/(N-1);
y=ymin:hy:ymax;

h=hx;
% make 2d spatial grid
[X,Y]=ndgrid(x,y);


% define the initial displacement and velocity vs. z
% for parts a-c
z=exp(-5*(X.^2+Y.^2));
vz=zeros(length(x),length(y));
% for part d
% vz=exp(-5*(X.^2+Y.^2));
% z=zeros(length(x),length(y));
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
% pause;

% Choose how long to run and when to plot
% tfinal=input(' Enter tfinal - ');
tfinal=100;
skip=input(' Enter # of steps to skip between plots (faster) - ');
nsteps=tfinal/tau;

% here is the loop that steps the solution along

% define an array for z(0,0,t)
zZeros=zeros(1,floor(nsteps));
close all;
figure  % open a new frame for the animation
xlabel('x');
ylabel('y');
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
      title(['Staggered Leapfrog Wave: Time= ' num2str(time)])
      caxis([-.25 .25])
      pause(.1)
   end
end


% plot z(0,0,t)
figure
times=tau:tau:tfinal;
plot(times,zZeros);
title('Z(0,0,t)');
xlabel('time');
ylabel('z(t)');

% Take the FFT of this trace and plot in a new figure window.
zN=floor((tfinal-tau)/tau);
freqs2 = tau*(0:(zN/2)); % Includes aliasing
freqs1 = freqs2(1:floor(zN/2+1)); 
P2 = abs(fft(zZeros)/zN);
P1 = P2(1:floor(zN/2+1));
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(freqs1,P1);
axis([0 2 0 max(P1(1:floor(2/tau)))*1.2]); % Adjust the axes so that you can see the range of 0-2 Hz.
% (Note that this is Hz, not radians/second, so we want things plotted in
% terms of f, not ?.)
title('Fourier Spectrum')
xlabel('f (Hz)')
ylabel('Amplitude')
% The spikes on these plots represent the frequencies of
% the normal modes of the wave equation in this geometry. 
% After making the plot of the FFT, have your script write the lowest 5 peaks
% visible in the plot and the (n, m) pairs that correspond to those peaks.
[peaks,locations] = findpeaks(P1,freqs1);
maxnm=25;
nn=ones(maxnm,1)*(1:1:maxnm);
mm=transpose((1:1:maxnm))*ones(1,maxnm);
freqpeaks=c/(4*L)*sqrt(nn.^2 + mm.^2);
for n1=1:5
    fprintf('peak #%g approximately at %f Hz\n',n1,locations(n1));
    indeces = abs(locations(n1)-freqpeaks)./freqpeaks < .015;
    % If there are degenerate peaks (peaks that are shared by more than one
    % combination of n and m) indicate all of the combinations that produce
    % that peak.
    for n2=1:1:maxnm
        for m2=1:1:maxnm
            if (indeces(n2,m2))
                fprintf('\tThis corresponds with a theoretical frequency of %f at n=%g and m=%g\n',freqpeaks(n2,m2),n2,m2);
            end;
        end;
    end;
end;
% The initial displacement and velocity are responsible for which modes are
% excited in this scenario. The reason why not all modes are used is
% because they're not included in this particular initial condition.
pause;


















% START OF PROBLEM 2
clear;close all;

% Set the number of grid points and build a cell-center grid with ghost
% points
scaler=1;
N=input(' Enter N, cell number - ');
L=1*scaler;
h=L/N;
x=-.5*h:h:L+.5*h;
x=x';  % Turn x into a column vector.

%  Load the diffusion coefficient array (make it a column vector)
D=2; % (just 2 for now--we'll change it later)

% Set the initial temperature distribution
% T=sin(pi*x/L);
T=ones(N+2)*10;

% Find the maximum of T for setting plot limits
Tmax=12;Tmin=0;

% Choose the time step tau.
% The max tau for explicit stability is a reasonable choice
% fprintf(' Maximum explicit time step: %g \n',h^2/max(D));
% tau = input(' Enter the time step - ');
tau=h^2/max(D);

% Create the matrices A and B by loading them with zeros
A=zeros(N+2);
B=zeros(N+2);

% load A and B at interior points
const = 2*h^2 / tau;
for j=2:N+1
   A(j,j-1)= -(1/(2*h^2));
   A(j,j)  = (1/(tau*D)+1/(2*x(j)*h)+1/(h.^2));
   A(j,j+1)= -(1/(2*h^2)+1/(2*x(j)*h));

   B(j,j-1)= (1/(2*h^2));
   B(j,j)  = (1/(tau*D)-1/(2*x(j)*h)-1/(h.^2));
   B(j,j+1)= (1/(2*h^2)+1/(2*x(j)*h));
end

% load the boundary conditions into A and B
A(1,1)=1/h; A(1,2)=-1/h; B(1,1)=0.;      % T'(0)=0
A(N+2,N+1)=0.5; A(N+2,N+2)=0.5; B(N+2,N+2)=0; % T(L)=0

% Set the number of time steps to take.
% tfinal=input(' Enter the total run time - ');
tfinal=.5;
nsteps=tfinal/tau;

% Choose how many iterations to skip between plots
% skip = input(' Number of iterations to skip between plots - ');
skip=1;
% This is the time advance loop.
for mtime=1:ceil(nsteps)
    % define the time
    t=mtime*tau;
    % find the right-hand side for the solve at interior points
    r=B*T;

    % apply the boundary conditions
    r(1)=0;   % T'(0)=0
    r(N+2)=0; % T(L)=0

    % do the linear solve to update T
    T=A\r;

    % Make a plot of T every once in a while.
    if(rem(mtime,skip) == 0)
        plot(x/scaler,T);
        axis([0 L/scaler Tmin Tmax])
        title(sprintf('Approximate Temperature: Time= %f',t))
        ylabel('Temp (C)');
        xlabel('rho');
        pause(.001)
    end
end
pause;
clear; close all;