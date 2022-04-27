load rf1.mat;
myData = rf1;
[nr,ns] = size(myData);
[x,y] = meshgrid(1:1:nr,1:1:ns);

Drr = myData(1:8:nr,1:4:ns);
[nrr,nss] = size(Drr);
[X,Y] = meshgrid(1:1:nrr,1:1:nss);

% x and y coordinates
xData = reshape(x,nr*ns,1) ;
yData = reshape(y,nr*ns,1) ;

% c data for the (x,y) points
cData = reshape(myData,nr*ns,1) ;

% x and y limits
xLim = [1 nrr] ;
yLim = [1 nss] ;

% % regular scatter plot
% figure
% scatter(xData,yData,20,cData,'marker','*')
% title('Regular scatter plot')
% axis equal ; box on
% colorbar ;

% hexScatter plot with rHex = 1.0 without filling empty hexagons
rHex = 1.0 ;
points = hexScatterm(xData,yData,cData,xLim,yLim,rHex) ;
pn = reshape(points,nr,ns) ;

