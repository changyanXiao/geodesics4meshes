function Trig=MakeTrig(Pts,Facets,Edges)

Trig.Coords = Pts;
Trig.Facets = Facets;
Trig.Edges  = Edges;
Trig.FixIds = [];

Trig.NVerts = size(Trig.Coords,1);
Trig.NFacets= size(Trig.Facets,1);