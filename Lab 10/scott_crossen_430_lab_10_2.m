clear; close all;
oldx=0;
newx=.5;
while(abs((newx-oldx)/(newx))>.0001)
    oldx=newx;
    newx=exp(-newx);
end;
fprintf('x= %f\n',newx)

oldx=0;
newx=.5;
for i=1:1:100
    oldx=newx;
    newx=-log(newx);
end;
fprintf('x= %f after %g iterations\n',newx,i);
