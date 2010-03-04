function Ls = computeTrigEdgeLengths(Trig)

N = size(Trig.Edges,2);
Ls = zeros(N,1);

for i = 1:N
   i1 = Trig.Edges(1,i); 
   i2 = Trig.Edges(2,i);
   
   c1 = Trig.Coords(i1,:);
   c2 = Trig.Coords(i2,:);
   
   Ls(i) = norm(c1-c2);
end