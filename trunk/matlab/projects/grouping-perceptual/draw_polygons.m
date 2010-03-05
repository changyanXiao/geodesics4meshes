function sk = draw_polygons(n,r,point_list)

if not(iscell(point_list))
    point_list = {point_list};
end

sk = zeros(n);
for i=1:length(point_list)
    pl = point_list{i};
    for k=2:length(pl)
        sk = draw_line(sk,pl(1,k-1),pl(2,k-1),pl(1,k),pl(2,k),r);
    end
end


function sk = draw_line(sk,x1,y1,x2,y2,r)


n = size(sk,1);
[Y,X] = meshgrid(1:n,1:n);
q = 80;
t = linspace(0,1,q);
x = x1*t+x2*(1-t); y = y1*t+y2*(1-t);
if r==0
    x = round( x ); y = round( y );
    sk( x+(y-1)*n ) = 1;
else
    for k=1:q
        I = find((X-x(k)).^2 + (Y-y(k)).^2 <= r^2 );
        sk(I) = 1;
    end
end