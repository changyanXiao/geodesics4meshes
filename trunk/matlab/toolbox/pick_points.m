%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick_points.m : Points = pick_points(Image, nb_points, msg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = pick_points(n)
    
    P = ones(2,n);
    
    for nbp = 1:n
        start_point = ginput(1);
        plot(start_point(1),start_point(2),'rs');
        P(1,nbp)=start_point(2);
        P(2,nbp)=start_point(1);
    end
    P = round(P);
    
end
