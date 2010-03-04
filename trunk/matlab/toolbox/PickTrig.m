% Pick a vertex of an edge of a triangulation overlaid on an image. 
function PickTrig(Trig,Proj,Img,FacetsP,EdgesP)

if(nargin<3)
    Img=[];
end
if(nargin<4)
    FacetsP=false;
end
if(nargin<5)
    EdgesP=true;
end

% Open a new figure. 
figure;

% Compute vertex projections
[Us,Vs]=ProjectTrig(Trig,Proj);
Us=round(Us);
Vs=round(Vs);
% Overlay edges on image and create a raster array. 
PickIm=FillIm(Trig,Us,Vs,Proj,Img,FacetsP,EdgesP);
% Capture mouse clicks. 
ClickIm(Trig,PickIm,Us,Vs);

end

%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

function ClickIm(Trig,PickIm,Us,Vs)
% Displays an image with vertices and edges superimposed on top of it.
% Implements the cursor data update function to show additional imformation
% on clicked vertices and edges.

% Get handle to current figure. 
h=get(0,'CurrentFigure');

daspect([1,1,1]);
datacursormode on;
hold on;
dcm_obj = datacursormode(h);
set(dcm_obj,'UpdateFcn',@CursorDataUpdateFcn,'SnapToDataVertex','off')

v1=0;
v2=0;
v3=0;

    function  CursorText = CursorDataUpdateFcn(empt, eventdata)
        
        Pos    = get(eventdata,'Position');
        %ImgObj = get(eventdata,'Target');
        CursorText = {''};
        %disp(size(PickIm));
        %disp(Pos);
        x=round(Pos(1));
        y=round(Pos(2));
        if((0<x)&&(0<y)&&(x<=size(PickIm,1))&&(y<=size(PickIm,2)))
            Id=PickIm(x,y);
            %fprintf('CursorDataUpdateFcn: %d %d -> %d\n',x,y,Id);
        else
            fprintf('Clicked outside of image.\n');
            Id=0;
        end
        if(Id>0)
            % Redraw in blue any edge that was overlaid in red.
            if((v1>0)&&(v2>0))
                plot([Us(v1),Us(v2)],[Vs(v1),Vs(v2)],'b');
                if(v3>0)
                    plot([Us(v1),Us(v3)],[Vs(v1),Vs(v3)],'b');
                    plot([Us(v2),Us(v3)],[Vs(v2),Vs(v3)],'b');
                    v3=0;
                end
                v1=0;
                v2=0;
            end
            if(1==mod(Id,3))
                % An edge has been selected.
                EId=fix(Id/3);
                v1=Trig.Edges(1,EId);
                v2=Trig.Edges(2,EId);
                % Draw selected edge in red. 
                plot([Us(v1),Us(v2)],[Vs(v1),Vs(v2)],'r');
                CursorText  = {  ['V1:  ',num2str(v1)],['V2:  ',num2str(v2)],['Id: ',num2str(EId)] };
                
            elseif(2==mod(Id,3))
                % A vertex has been selected.
                VId=fix(Id/3);
                x=Trig.Coords(VId,1);
                y=Trig.Coords(VId,2);
                z=Trig.Coords(VId,3);
                CursorText  = {  ['X:  ',num2str(x)],['Y:  ',num2str(y)],['Z:  ',num2str(z)],['Id: ',num2str(VId)] };
            else
                % A facet has been selected.
                FId=fix(Id/3);
                v1=Trig.Facets(FId,1);
                v2=Trig.Facets(FId,2);
                v3=Trig.Facets(FId,3);
                plot([Us(v1),Us(v2)],[Vs(v1),Vs(v2)],'r');
                plot([Us(v1),Us(v3)],[Vs(v1),Vs(v3)],'r');
                plot([Us(v2),Us(v3)],[Vs(v2),Vs(v3)],'r');
                CursorText  = {  ['V1:  ',num2str(v1)],['V2:  ',num2str(v2)],['V3:  ',num2str(v3)],['Id: ',num2str(FId)] };
                
            end
        end
        %disp(CursorText);
    end
end

%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

function [PickIm,DispIm]=FillIm(Trig,Us,Vs,Proj,Img,FacetsP,EdgesP)

if(isempty(Img))
    % No Img is specified
    xdim=max(Us);
    ydim=max(Vs);
    Img=zeros(ydim,xdim);
else
    % Img is specified
    ydim=size(Img,1);
    xdim=size(Img,2);
end
% PickIm and DispIm are raster arrays and therefore transposed.
PickIm=zeros(xdim,ydim);
% Fill PickIm with FacetIds if the facets are supposed to be sensitive.
if(FacetsP)
    PickIm=3*FillTrig(Trig,Proj,PickIm,true);
end
if(nargout>1)
    DispIm=zeros(xdim,ydim);
else
    DispIm=[];
end

VertVal=1.0;
EdgeVal=1.0;

clf;
imshow(Img);
hold on;

% Draw scan lines between edges when both vertices are in the image.
for EId = 1:length(Trig.Edges)
    i1 = Trig.Edges(1,EId);
    i2 = Trig.Edges(2,EId);
    x1= Us(i1);
    y1= Vs(i1);
    x2= Us(i2);
    y2= Vs(i2);
    %fprintf('%f %f %f %f:\n',x1,y1,x2,y2);
    if((0<x1)&&(0<y1)&&(x1<=xdim)&&(y1<=ydim)&&(0<x2)&&(0<y2)&&(x2<=xdim)&&(y2<=ydim))
        plot([x1 x2],[y1 y2],'b');
        % Scan line from (x1,y1) to (x2,y2)
        if(EdgesP)
            [xs,ys]=ScanLine(x1,y1,x2,y2);
            for j=2:length(xs)-1
                x=xs(j);
                y=ys(j);
                v=3*EId+1;
                PickIm(x,y)=v;
                if(x<xdim)
                    PickIm(x+1,y)=v;
                end
                if(y<xdim)
                    PickIm(x,y+1)=v;
                end
                if(1<x)
                    PickIm(x-1,y)=v;
                end
                if(1<y)
                    PickIm(x,y-1)=v;
                end
                if(~isempty(DispIm))
                    DispIm(x,y)=EdgeVal;
                end
            end
        end
    end
end
% For each vertex, draw a cross that are five pixel wide in both horizontal 
% and vertical directions.
for VId = 1:length(Us)
    x = Us(VId);
    y = Vs(VId);
    v = 3*VId+2;
    if((0<x)&&(0<y)&&(x<=xdim)&&(y<=ydim))
        plot(x,y,'.r');
        PickIm(x,y)=v;
        if(~isempty(DispIm))
            DispIm(x,y)=VertVal;
        end
        if(x<xdim)
            PickIm(x+1,y)=v;
            if(~isempty(DispIm))
                DispIm(x+1,y)=VertVal;
            end
            if(x<xdim-1)
                PickIm(x+2,y)=2*VId;
                if(~isempty(DispIm))
                    DispIm(x+2,y)=VertVal;
                end
            end
         end
        if(y<ydim)
            PickIm(x,y+1)=v;
            if(~isempty(DispIm))
                DispIm(x,y+1)=VertVal;
            end
            if(y<ydim-1)
                PickIm(x,y+2)=v;
                if(~isempty(DispIm))
                    DispIm(x,y+2)=VertVal;
                end
            end
        end
        if(1<x)
            PickIm(x-1,y)=v;
            if(~isempty(DispIm))
                DispIm(x-1,y)=VertVal;
            end
            if(2<x)
                PickIm(x-2,y)=v;
                if(~isempty(DispIm))
                    DispIm(x-2,y)=VertVal;
                end
            end
        end
        if(1<y)
            PickIm(x,y-1)=v;
            if(~isempty(DispIm))
                DispIm(x,y-1)=VertVal;
            end
            if(1<y)
                PickIm(x,y-2)=v;
                if(~isempty(DispIm))
                    DispIm(x,y-2)=VertVal;
                end
            end
        end
    end
end
hold off;
end

% function [Us,Vs]=ProjectTrig(Trig,Proj)
% 
% assert((3==size(Proj,1))&&(4==size(Proj,2)));
% XYZs=[Trig.Coords';ones(1,size(Trig.Coords,1))];
% UVWs=Proj*XYZs;
% Us=round(UVWs(1,:)./UVWs(3,:));
% Vs=round(UVWs(2,:)./UVWs(3,:));
% 
% end