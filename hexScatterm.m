%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE REDUCES THE DATA RESOLUTION BY PREPARING A HEXAGONAL 2D GRID %  
% AND REPLACING ALL THE DATA POINTS INSIDE A HEXAGONE WITH A SINGLE      %
% VALUE. IT THEN GIVES OUT A SCATTER PLOT WITH THOSE HEXAGONES.          %
% PREPARED BY: SALMAN MASHAYEKH                                          %
% DATE: DECEMBER 2012                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hexMmbr = hexScatterm(xData,yData,cData,xLim,yLim,rHex)

% make sure the inputs are vertical
if ~isvector(xData) || ~isvector(yData) || ~isvector(cData) || ...
        length(xData)~=length(yData) || length(xData)~=length(cData) || length(yData)~=length(cData)
    fprintf('Error! xData, ydata, and cData must be vectors of the same length!\n');
    return
end
xData = reshape(xData,length(xData),1) ;
yData = reshape(yData,length(yData),1) ;
cData = reshape(cData,length(cData),1) ;

% make the grid for the center of hexagons
xrow1 = xLim(1)+sqrt(3)/2*rHex:sqrt(3)*rHex:xLim(2)+sqrt(3)/2*rHex ;
xrow2 = xLim(1):sqrt(3)*rHex:xLim(2) ;
if length(xrow1) > length(xrow2)
    xrow1 = xrow1(1:end-1) ;
elseif length(xrow2) > length(xrow1)
    xrow2 = xrow2(1:end-1) ;
end
yrow1 = yLim(1):3/2*rHex:yLim(2) ;
% ------
xgrid = repmat([xrow2;xrow1],floor(length(yrow1)/2),1) ;
if mod(length(yrow1),2) == 1
    xgrid = [xrow1; xgrid] ;
end
ygrid = repmat(yrow1(1:size(xgrid,1))',1,length(xrow1)) ;
% ------
xHex = xgrid(:) ;
yHex = ygrid(:) ;

% determine which hexagon each (x,y) data belongs to
hexMmbr = zeros(size(xData)) ;
for n = 1:length(xData)
    distVect = (xHex-xData(n)).^2+(yHex-yData(n)).^2 ;
    [~,hexMmbr(n)] = min(distVect) ;
end
