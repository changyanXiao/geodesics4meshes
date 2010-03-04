% Project vertex by multiplying coordinates by projection matrix and
% adding one to convert to matlab array coordinates. 
function [u,v]=ProjectVert(Trig,Proj,VId,RoundP)

if(nargin<4)
    RoundP=false;
end

Pt=Trig.Coords(VId,:);
[u,v]=ProjectPoint(Proj,Pt',RoundP);

