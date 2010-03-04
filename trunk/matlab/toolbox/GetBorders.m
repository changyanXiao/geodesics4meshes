% if OrderP=true, ensure that the edges are ordered so that the second of
% one is the first vertex of the next one.
% If ReverseP=true, reverse the order of the vertices in all edges. 
function [Borders,Edges]=GetBorders(Trig,OrderP,ReverseP)

if(nargin<3)
    ReverseP=false;
    if(nargin<2)
        OrderP=true;
    end
end

nv=size(Trig.Coords,1);
ne=0;
nb=0;

Borders=zeros(2,5*nv);
if(nargout>1)
    Edges=zeros(2,5*nv);
end
Count=zeros(nv,nv);

if(OrderP)
    % Make an ordered list of edges that form the order.
    for i=1:size(Trig.Facets,1)
        j1=Trig.Facets(i,1);
        j2=Trig.Facets(i,2);
        j3=Trig.Facets(i,3);
        Count(j1,j2)=Count(j1,j2)+1;
        Count(j2,j1)=Count(j2,j1)+1;
        Count(j1,j3)=Count(j1,j3)+1;
        Count(j3,j1)=Count(j3,j1)+1;
        Count(j2,j3)=Count(j2,j3)+1;
        Count(j3,j2)=Count(j3,j2)+1;
    end
    % Find all edges if a second output argument is specified. 
    if(nargout>1)
        for j=1:nv
            for i=1:j-1
                if(Count(i,j)>0)
                    ne=ne+1;
                    Edges(1,ne)=i;
                    Edges(2,ne)=j;
                end
            end
        end
        Edges=Edges(:,1:ne);
    end
    % Find border edges
    for j=1:nv
        for i=1:j-1
            % Find the first edge and store the index of its first vertex.
            if(Count(i,j)==1)
                Count(i,j)=0;
                Count(j,i)=0;
                if(ReverseP)
                    StrtI=j;
                    NextI=i;
                else
                    StrtI=i;
                    NextI=j;
                end
                nb=nb+1;
                Borders(1,nb)=StrtI;
                Borders(2,nb)=NextI;
                % Follow the boundary until reaching the first vertex again.
                while(NextI~=StrtI)
                    for k=1:nv
                        if(Count(NextI,k)==1)
                            Count(NextI,k)=0;
                            Count(k,NextI)=0;
                            nb=nb+1;
                            Borders(1,nb)=NextI;
                            Borders(2,nb)=k;
                            NextI=k;
                            break;
                        end
                    end
                end
                % Adjust array size and return.
                Borders=Borders(:,1:nb);
                return;
            end
        end
    end
else
    % Find the border edges in random order and optionally also find all
    % the other edges.
    for i=1:size(Trig.Facets,1)
        j1=Trig.Facets(i,1);
        j2=Trig.Facets(i,2);
        j3=Trig.Facets(i,3);
        if(j1<j2)
            Count(j1,j2)=Count(j1,j2)+1;
        else
            Count(j2,j1)=Count(j2,j1)+1;
        end
        if(j1<j3)
            Count(j1,j3)=Count(j1,j3)+1;
        else
            Count(j3,j1)=Count(j3,j1)+1;
        end
        if(j2<j3)
            Count(j2,j3)=Count(j2,j3)+1;
        else
            Count(j3,j2)=Count(j3,j2)+1;
        end
    end
    
    if(nargout>1)
        for j=1:nv
            for i=1:j-1
                if(Count(i,j)>0)
                    ne=ne+1;
                    Edges(1,ne)=i;
                    Edges(2,ne)=j;
                    if(Count(i,j)==1)
                        nb=nb+1;
                        Borders(1,nb)=i;
                        Borders(2,nb)=j;
                    end
                end
            end
        end
    else
        for j=1:nv
            for i=1:j-1
                if(Count(i,j)==1)
                    nb=nb+1;
                    Borders(1,nb)=i;
                    Borders(2,nb)=j;
                end
            end
        end
    end
    Borders=Borders(:,1:nb);
    if(nargout>1)
        Edges=Edges(:,1:ne);
    end
end



    
    