function [M, nb_labels] = interactive_labeling(Mtx)
%construct label image mask

    bound = [];
    nb_labels = 1;
    n = size(Mtx,1);
    m = size(Mtx,2);
    M = zeros(n,m);
    max_label = 100;
    c = hsv(max_label);
    c_cur = c(1,:);
    figure;
    imagesc(Mtx);
    axis image; axis on;
    hold on;
    title('Click [left] for new point, [middle] to another polyline of same label, [right] to next label, [q] to stop labeling, [c] to close the polyline');
    while 1
        if not(isempty(bound))
            %hold on;
            h = plot(bound(2,:), bound(1,:), '.-', 'Color', c_cur);
            set(h, 'LineWidth', 2);
        end
        [x,y,b] = ginput(1);
        if b==3  % another class + curve
            if not(isempty(bound))
                M = draw_polygons(M, 2, bound, nb_labels);
                nb_labels = nb_labels + 1;
                if (nb_labels <= max_label)
                    c_cur = c(nb_labels,:);
                end
                bound = [];
            end
            continue;
        end  % another curve
        if b==2
            M = draw_polygons(M, 2, bound, nb_labels);
            bound = [];
            continue;
        end
        if b==1  % new point
            bound(:,end+1) = [y;x];
            continue;
        end
        if b==99  % close the polyline and go to the next polyline
            bound(:,end+1) = bound(:,1);
            h = plot(bound(2,end-1:end), bound(1,end-1:end), '.-', 'Color', c_cur);
            set(h, 'LineWidth', 2);
            M = draw_polygons(M, 2, bound, nb_labels);
            bound = [];
            continue;
        end
        if b==113  % quit the labeling process
            if not(isempty(bound))
                M = draw_polygons(M, 2, bound, nb_labels);
            end
            break;
        end
    end
    hold off;
end
