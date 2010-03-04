function DispNumberedTrig(Trig,vertexIndeciesP,faceIndeciesP,varargin)
figure
if(length(varargin)>0)
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
    
end

daspect([1,1,1]);
view(10,20);

hold on;
for i=1:length(Trig)
    Color=Colors(1+mod(i-1,Ncolrs));
    DrawTrig(Trig(i),Color,vertexIndeciesP,faceIndeciesP);
end
xlabel('X');
ylabel('Y');
zlabel('Z');
hold off;


%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

function DrawTrig(Trig,Color,vertexIndeciesP,faceIndeciesP)

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

% draw vertex indecies
if(vertexIndeciesP)
    for i=1:1:Trig.NVerts
        text(Trig.Coords(i,1),Trig.Coords(i,2),Trig.Coords(i,3),sprintf('%d',i))
    end
end

% draw vertex indecies
if(vertexIndeciesP)
    for i=1:1:Trig.NVerts
        text(Trig.Coords(i,1),Trig.Coords(i,2),Trig.Coords(i,3),sprintf('%d',i))
    end
end
% draw face indecies
if(faceIndeciesP)
    for i=1:1:Trig.NFacets
         i1 = Trig.Facets(i,1);i2 = Trig.Facets(i,2);i3 = Trig.Facets(i,3);   
         center(1) = sum(Trig.Coords(i1,1)+Trig.Coords(i2,1)+Trig.Coords(i3,1))/3;
         center(2) = sum(Trig.Coords(i1,2)+Trig.Coords(i2,2)+Trig.Coords(i3,2))/3;
         center(3) = sum(Trig.Coords(i1,3)+Trig.Coords(i2,3)+Trig.Coords(i3,3))/3;
         text(center(1),center(2),center(3),sprintf('%d',i))
    end
end
