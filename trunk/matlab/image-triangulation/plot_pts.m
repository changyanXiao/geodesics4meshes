function plot_pts(P,c,s)

if size(P,1)==2
    h = plot(P(2,:),P(1,:),c);
    set(h,'MarkerSize', s);
    return;
end
if size(P,1)==3
    h = plot3(P(3,:),P(2,:),P(1,:),c);
    set(h, 'MarkerSize', s);
    return;
end

end
