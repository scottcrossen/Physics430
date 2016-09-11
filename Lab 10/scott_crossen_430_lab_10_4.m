clear; close all;
oldx=0;
newx=.5;
omega=input('input an omega: ');
iteration=0;
% We can't solve the equation in 10.4 in one step because we don't take the
% series out that many parts.
while(abs((newx-oldx)/(newx))>.00001 && iteration <100) % reduce acceptable error to increase speed
    iteration=iteration+1;
    oldx=newx;
    newx=omega*exp(-newx)+(1-omega)*newx;
end;
fprintf('x= %f in %g iterations\n',newx,iteration)
