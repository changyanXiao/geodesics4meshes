function facets = getFacetsForVertex(Trig,vin)

% facets = facetId , localIndexOfVertex
% localIndexOfVertex = 1,2 or 3

facets = [];

firsts = Trig.Facets(:,1);
seconds = Trig.Facets(:,2);
thirds = Trig.Facets(:,3);

fins1 = find( firsts==vin);
fins2 = find( seconds==vin);
fins3 = find( thirds==vin);

for j = 1 : length(fins1)
    facets = [facets;[fins1(j),1]];
end

for j = 1 : length(fins2)
    facets = [facets;[fins2(j),2]];
end

for j = 1 : length(fins3)
    facets = [facets;[fins3(j),3]];
end



