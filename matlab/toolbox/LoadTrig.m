function Trig=LoadTrig(BaseName)

Trig.Coords = load(sprintf('%s.pts',BaseName));
Trig.Facets = load(sprintf('%s.tri',BaseName))+1;
%Trig.FixIds = LoadIfExists(sprintf('%s.fix',BaseName));

Trig.NVerts = size(Trig.Coords,1);
Trig.NFacets= size(Trig.Facets,1);

Trig.Edges  = GetEdges(Trig);
