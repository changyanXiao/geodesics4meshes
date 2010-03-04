% Project point by multiplying coordinates by projection matrix and
% adding one to convert to matlab array coordinates. 
function [u,v]=ProjectPoint(Proj,Pt,RoundP)

UVW=Proj*[Pt;1];
u=UVW(1)/UVW(3)+1;
v=UVW(2)/UVW(3)+1;

if((nargin>2)&&RoundP)
    u=round(u);
    v=round(v);
end
  