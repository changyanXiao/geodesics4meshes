function [eins] = getEdgesForVertex(Trig0,i)

edges = Trig0.Edges;

e1 = edges(1,:);
e2 = edges(2,:);

l1 = find(e1==i);
l2 = find(e2==i);

eins = [l1,l2];






