L=100;
tau=.01;
N=L/tau;
x=0:tau:L;
y=sin(x)+sin(3*x)+sin(5*x);
figure;
plot(x,y);
title('time domain')
freqs = tau*(0:(N/2));
P2 = abs(fft(y)/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(freqs*2*pi,P1);
axis([0 6 0 1]);
title('w (Hz)')
figure;
plot(freqs,P1);
axis([0 6/(2*pi) 0 1]);
title('f (Hz)')