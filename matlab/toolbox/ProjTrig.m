function ProjTrig(Trig,PrjMat,Img,RasterP,Colors,ClearP)


if(nargin<6)
    ClearP=true;
    if(nargin<5)
        Colors=['r','g','b','c','y'];
        if(nargin<4)
            RasterP=false;
            if(nargin<3)
                Img=[];
            end
        end
    end
end
Ncolrs=length(Colors);

if(isempty(Img))
    if(ClearP)
        daspect([1 1 1]);
        clf;
    end
elseif(RasterP)
    DispRaster(Img,true);
else
    imshow(Img);
end

hold on;
for i=1:length(Trig)
    Color=Colors(1+mod(i-1,Ncolrs));
    DrawTrig(Trig(i),PrjMat,Color);
end
hold off;

%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

function DrawTrig(Trig,PrjMat,Color)

[Us,Vs]=ProjectTrig(Trig,PrjMat);

% Draw triangulation edges
if(isfield(Trig,'Edges'))
    for e=1:size(Trig.Edges,2)
        i=Trig.Edges(1,e);
        j=Trig.Edges(2,e);
        plot([Us(i),Us(j)],[Vs(i),Vs(j)],Color);
    end
elseif(isfield(Trig,'Neighbors'))
    for i=1:Trig.NVerts
        for j=i:Trig.NVerts
            if Trig.Neighbors(i,j)==1
                plot([Us(i),Us(j)],[Vs(i),Vs(j)],Color);
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
        plot(Us(j),Vs(j),'xb');
    end
end
