function [facetIds]=getFacetsForEdge(Trig,e)

%facetIds = facetId, localindexofstart, localindexofend
facetIds = [];

e1 = Trig.Edges(1,e);
e2 = Trig.Edges(2,e);

for i = 1 : Trig.NFacets
    
    if(Trig.Facets(i,1)==e1 && Trig.Facets(i,2)==e2)
       facetIds = [facetIds;[i,1,2]];
    end
    if(Trig.Facets(i,1)==e2 && Trig.Facets(i,2)==e1)
       facetIds = [facetIds;[i,2,1]];
    end
    if(Trig.Facets(i,1)==e1 && Trig.Facets(i,3)==e2)
        facetIds = [facetIds;[i,1,3]];
    end
    if(Trig.Facets(i,1)==e2 && Trig.Facets(i,3)==e1)
        facetIds = [facetIds;[i,3,1]];
    end
    if(Trig.Facets(i,2)==e1 && Trig.Facets(i,3)==e2)
        facetIds = [facetIds;[i,2,3]];
    end
    if(Trig.Facets(i,2)==e2 && Trig.Facets(i,3)==e1)
        facetIds = [facetIds;[i,3,2]];
    end
    
end