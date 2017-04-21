function [ newXvect, newYvect ] = AK_rotateTranslate2D( Xvect, Yvect, rotationAngle, Xcenter, Ycenter )
%AK_rotateTranslate2D takes xy-coordinates about the origin for some shape
%(i.e. patch), rotates the shape, and translates it to a new center position
%   INPUT
%       Xvect: vector of x coordinates for shape centered on origin
%       Yvect: vector of y coordinates for shape centered on origin
%       rotationAngle: angle of rotation
%       Xcenter: new center position for shape (x)
%       Ycenter: new center position for shape (y)
%   OUTPUT
%       newXvect: vector of new x coordinates for shape
%       newYvect: vector of new y coordinates for shape


% rotation matrix
rotMat = [cos(rotationAngle) sin(rotationAngle) 0; -sin(rotationAngle) cos(rotationAngle) 0; 0 0 1];

% translation matrix
transMat = [1 0 0; 0 1 0; Xcenter Ycenter 1];

% preallocate
newXvect = nan(length(Xvect),1);
newYvect = nan(length(Yvect),1);
for iCoord = 1:length(Xvect)
    % operate on temporary coodinates in a homogeneous coordinate space
    clear temp
    % rotate
    temp = [Xvect(iCoord) Yvect(iCoord) 1]*rotMat;
    %translate to new position
    temp = [temp(1) temp(2) 1]*transMat;
    % assignment
    newXvect(iCoord) = temp(1);
    newYvect(iCoord) = temp(2);
end

