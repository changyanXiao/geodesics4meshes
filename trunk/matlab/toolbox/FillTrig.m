function [Raster,ValidP]=FillTrig(Trig,Proj,Raster,FacetsP)

if(nargin<4)
    FacetsP=false;
end

global Array;
Array=Raster;

if(FacetsP)
    % Fill array with individual facet numbers
    nf=size(Trig.Facets,1);
    % Compute vertex projections.
    [Us,Vs]=ProjectTrig(Trig,Proj);
    % Draw all edges into array.
    for EdgeId=1:size(Trig.Edges,2)
        n1=Trig.Edges(1,EdgeId);
        n2=Trig.Edges(2,EdgeId);
        DrawLine(Us(n1),Vs(n1),Us(n2),Vs(n2));
    end
    % Fill individual facets.
    us=zeros(nf,1);
    vs=zeros(nf,1);
    FIds=1:nf;
    for FacetId=FIds;
        n1=Trig.Facets(FacetId,1);
        n2=Trig.Facets(FacetId,2);
        n3=Trig.Facets(FacetId,3);
        us(FacetId)=round((Us(n1)+Us(n2)+Us(n3))/3.0);
        vs(FacetId)=round((Vs(n1)+Vs(n2)+Vs(n3))/3.0);
        
    end
    Raster=FillPoly(Array,FIds',us,vs);
    ValidP=true;
else
    % Same value for all facets.
    Pt=[mean(Trig.Coords(:,1));mean(Trig.Coords(:,2));mean(Trig.Coords(:,3))];
    [uc,vc]=ProjectPoint(Proj,Pt,true);
    if((0<uc)&&(0<vc)&&(uc<=size(Array,1))&&(vc<=size(Array,2)))
        if(Array(uc,vc))
            ValidP=false;
            return;
        end
    else
        ValidP=false;
        return;
    end
    % If the borders are already stored, use them. Otherwise recompute
    % them.
    if(isfield(Trig,'Borders'))
        Borders=Trig.Borders;
    else
        fprintf('Computing the triangulation borders.\n');
        Borders=GetBorders(Trig,true,false);
    end
    
    %clf;
    %hold on;
    for i=1:size(Borders,2)
        n1=Borders(1,i);
        n2=Borders(2,i);
        [u1,v1]=ProjectVert(Trig,Proj,n1,true);
        [u2,v2]=ProjectVert(Trig,Proj,n2,true); 
        %fprintf('%f %f %f %f\n',u1,v1,u2,v2);
        DrawLine(u1,v1,u2,v2,false);
        %plot([u1,u2],[v1,v2]);
    end
    %hold off;
    %pause;
    %DispRaster(Array);
    Raster=FillPoly(Array,1.0,uc,vc);
    ValidP=true;
end

function DrawLine(u1,v1,u2,v2,DbgP)


global Array;
%disp(size(Array));

xdim=size(Array,1);
ydim=size(Array,2);

[us,vs]=ScanLine(u1,v1,u2,v2);
%fprintf('%d %d: %d %d %d %d\n',xdim,ydim,u1,v1,u2,v2);
for j=1:length(us)
    u=us(j);
    v=vs(j);
    if((0<u)&&(0<v)&&(u<=xdim)&&(v<=ydim))
        %fprintf('.');
        Array(us(j),vs(j))=1;
    elseif(DbgP)
        fprintf('(%d,%d) out of bounds (%d,%d)\n',us(j),vs(j),xdim,ydim);
        return;
    end
end

%DispRaster(Array);
%pause;

