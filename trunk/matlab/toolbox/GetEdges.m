function Edges=GetEdges(Trig)
    
Neighbors=NeighborhoodFromNodes(Trig.Facets,Trig.NVerts);
Edges=[];
for j=1:Trig.NVerts
    for i=1:j-1
        if(Neighbors(i,j))
            Edges=[Edges,[i;j]];
        end
    end
end

function Neighborhood = NeighborhoodFromNodes(Nodes,nPts)

Neighborhood = zeros(nPts,nPts);
for i=1:size(Nodes,1)
    Neighborhood(Nodes(i,1),Nodes(i,2)) = 1;
    Neighborhood(Nodes(i,1),Nodes(i,3)) = 1;
    Neighborhood(Nodes(i,2),Nodes(i,1)) = 1;
    Neighborhood(Nodes(i,2),Nodes(i,3)) = 1;
    Neighborhood(Nodes(i,3),Nodes(i,1)) = 1;
    Neighborhood(Nodes(i,3),Nodes(i,2)) = 1;
end

