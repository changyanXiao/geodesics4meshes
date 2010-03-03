function plot_dependance_graph2D(point, dUx, dUy)



%[Nx,Ny] = size(dUx);
if dUx(point) > 0
    npoint1 = point-[1 0];
else if dUx(point) < 0
        npoint1 = point+[1 0];
    else npoint1= point;
    end;
end;
if dUy(point) > 0
    npoint2 = point-[0 1];
else if dUy(point) < 0
        npoint2 = point+[0 1];
    else npoint2 = point;
    end;
end;

if (sum(npoint1(:)-point(:)) ~=0)
    plot([point(2) npoint1(2)], [point(1) npoint1(1)], '-');
    plot_dependance_graph2D(npoint1, dUx, dUy);
end;
%if (sum(npoint2(:)-point(:)) ~=0)
%    plot([point(2) npoint2(2)], [point(1) npoint2(1)],'-');
%    plot_dependance_graph2D(npoint2, dUx, dUy);
%end;
