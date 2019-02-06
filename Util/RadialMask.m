function outputMask = RadialMask(numLines,X,Y,q)
%
% Compute radial mask
%
%
% Input:
%   numLines: number of radial lines through the origin
%   X: number of rows
%   Y: number of columns
%

imSize = max(X,Y);
m = ceil(sqrt(2)*imSize);
aux = zeros(m,m); ima = aux;
aux(round(m/2+1),:) = 1;

gamma = (sqrt(5)+1)/2;


angle = 180/gamma;

%angles = 0:angle:180-angle;
angles = angle.*[1+q*numLines:numLines+q*numLines];

for a = 1:length(angles)
    ang = angles(a);
    a = imrotate(aux,ang,'crop');
    ima = ima + a;
end
ima = ima(round(m/2+1) - imSize/2:round(m/2+1) + imSize/2-1,...
          round(m/2+1) - imSize/2:round(m/2+1) + imSize/2-1);

y = (ima > 0);

outputMask = double(y(imSize/2 +1 - X/2:imSize/2 + X/2,imSize/2 +1 - Y/2:imSize/2 + Y/2));
outputMask(round(X/2),round(Y/2)) = 1;

