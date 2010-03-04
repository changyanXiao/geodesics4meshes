% Project vertices by multiplying coordinates by projection matrix and
% adding one to convert to matlab array coordinates. 
function [Us,Vs]=ProjectTrig(Trig,Proj)

if((3==size(Proj,1))&&(4==size(Proj,2)))
    scale=1;
    XYZs=[Trig.Coords';ones(1,size(Trig.Coords,1))];
    UVWs=Proj*XYZs;
    Us=scale * UVWs(1,:)./UVWs(3,:);
    Vs=scale * UVWs(2,:)./UVWs(3,:);
elseif ((1==size(Proj,1))&&(1==size(Proj,2)))
    scale=Proj;
    Us=scale* Trig.Coords(:,1)./Trig.Coords(:,3);
    Vs=scale* Trig.Coords(:,2)./Trig.Coords(:,3);
else
    error('Proj should be either a 3x4 matrix or a scalar');
end

% Convert to matlab raster coordinates by adding one. 
Us=Us+1;
Vs=Vs+1;