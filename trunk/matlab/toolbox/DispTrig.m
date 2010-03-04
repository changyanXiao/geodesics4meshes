function DispTrig(Trig,varargin)

if(~isempty(varargin))
    Colors=varargin{1};
    if(length(varargin)>1)
        ClearP=varargin{2};
    else
        ClearP=true;    
    end
else
    ClearP=true;
    Colors=['r','g','b','c','y'];
end
Ncolrs=length(Colors);

if(ClearP)
    clf;
    daspect([1,1,1]);
end

hold on;
for i=1:length(Trig)
    Color=Colors(1+mod(i-1,Ncolrs));
    DrawTrig(Trig(i),Color);
end
hold off;

%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

function DrawTrig(Trig,Color)

% Draw triangulation edges
if(isfield(Trig,'Edges'))
    for e=1:size(Trig.Edges,2)
        i=Trig.Edges(1,e);
        j=Trig.Edges(2,e);
        plot3([Trig.Coords(i,1),Trig.Coords(j,1)],[Trig.Coords(i,2),Trig.Coords(j,2)],[Trig.Coords(i,3),Trig.Coords(j,3)],Color);
    end
elseif(isfield(Trig,'Neighbors'))
    for i=1:Trig.NVerts
        for j=i:Trig.NVerts
            if Trig.Neighbors(i,j)==1
                plot3([Trig.Coords(i,1),Trig.Coords(j,1)],[Trig.Coords(i,2),Trig.Coords(j,2)],[Trig.Coords(i,3),Trig.Coords(j,3)],Color);
            end
        end
    end
else
    warning('Cannot draw the edges of the triangulation');
end

% Draw cross at fixed vertices
if(isfield(Trig,'FixIds'))
    for i=1:length(Trig.FixIds)
        j=Trig.FixIds(i)+1;
        plot3(Trig.Coords(j,1),Trig.Coords(j,2),Trig.Coords(j,3),'xb');
    end
end
