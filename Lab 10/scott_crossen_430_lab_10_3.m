clear; close all;
oldx=0;
newx=.5;
omega=input('input an omega: ');
iteration=0;
while(abs((newx-oldx)/(newx))>.00001 && iteration <100)
    iteration=iteration+1;
    oldx=newx;
    newx=omega*exp(-newx)+(1-omega)*newx;
end;
fprintf('x= %f in %g iterations\n',newx,iteration)
