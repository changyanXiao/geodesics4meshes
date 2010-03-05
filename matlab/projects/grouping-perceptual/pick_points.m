%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick_points.m : Points = pick_points(Image, nb_points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = pick_points(M, n)
    
    P = ones(2,n);
    figure;
    imageplot(M);
    axis image; axis off;
    hold on;
    for nbp = 1:n
        disp('Pick a point.');
        start_point = ginput(1);
        plot(start_point(1),start_point(2),'rs');
        P(1,nbp)=start_point(2);
        P(2,nbp)=start_point(1);
    end
    hold off;
    P = round(P);
    
end
