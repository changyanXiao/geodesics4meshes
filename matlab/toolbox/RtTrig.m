% Rotate all vertex coordinates. 
% - If two or three arguments are provided, the second is assumed to be the
% rotations translation matrix and the third a flag that indicates whether
% or not it should be inverted. 
% - If six arguments are provided, they are assumed to be the rotation and
% translation values. 
function Trig=RtTrig(Trig,rx,ry,rz,tx,ty,tz)

if(2==nargin)
    RtMat=rx;
    assert(((3==size(RtMat,1))||(4==size(RtMat,1)))&&(4==size(RtMat,2)));
elseif(3==nargin)
    RtMat=rx;
    assert(((3==size(RtMat,1))||(4==size(RtMat,1)))&&(4==size(RtMat,2)));
    % Invert the RtMatrix
    if(ry)
        if(3==size(RtMat,1))
            RtMat=inv([RtMat;[0,0,0,1]]);
        else
            RtMat=inv(RtMat);
        end
    end
elseif (6==nargin)
    % Compute the RtMatrix from the three rotations and translations.
    RtMat=Rot4(rx,ry,rz,tx,ty,tz);
else
    error('RtTrig expects either 2, 3, or 6 elements');
end

RtPts=[Trig.Coords,ones(Trig.NVerts,1)];
RtPts=RtMat*RtPts';
RtPts=RtPts';
Trig.Coords=RtPts(:,1:3);