function SaveTrig(BaseName,Trig)

Pts=Trig.Coords;
save(sprintf('%s.pts',BaseName),'Pts','-ascii');
Tri=Trig.Facets-1;
file=fopen(sprintf('%s.tri',BaseName),'w');
for i=1:Trig.NFacets
    fprintf(file,'%d %d %d\n',Tri(i,1),Tri(i,2),Tri(i,3));
end
fclose(file);

if(isfield(Trig,'FixIds'))
    file=fopen(sprintf('%s.fix',BaseName),'w');
    for i=1:length(Trig.FixIds)
        fprintf(file,'%d\n',Trig.FixIds(i));
    end
    fclose(file);
end