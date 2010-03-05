%--------------------------------------------------------------------------
clear all; close all; clc;
%--------------------------------------------------------------------------
library_directory = '../../fmlib/';
addpath(library_directory);
data_directory   = '../../Data/';
addpath(data_directory);
matlab_directory   = '../../matlab/';
addpath(matlab_directory);
%--------------------------------------------------------------------------
DisplaySourcesColor=[1 0 0];
DisplayKeyPointsColor=[1 1 0];
DisplayMarkerSize=4;
%--------------------------------------------------------------------------
im = imread('Cavern.png');
Nx = 200; Ny = 200;
Ob = zeros(Nx+2, Ny+2);
Ob = logical(Ob);
Ob(2:Nx+1, 2:Ny+1) = (im < 100);
Obb = Ob(2:Nx+1, 2:Ny+1);
Pold = ones(Nx,Ny).*(1-Obb); clear Obb;
%--------------------------------------------------------------------------
w = 0.001;
h = [1;1];
source_points = [32;188];
destination_points = [186;18];
sz_source = size(source_points); nb_sources = sz_source(2);
sz_dest   = size(destination_points); nb_dest = sz_dest(2);
Traffic = ones(nb_sources, nb_dest);
Nb_iter = 400;
Energy = ones(1,Nb_iter-1);
N = 1;
Nvar = 0;
K = 5;
rho = 1/K;
%--------------------------------------------------------------------------
fig = figure(1);
%--------------------------------------------------------------------------
while(N < Nb_iter)
    %----------------------------------------------------------------------
    P = (1.0/12.0)*Pold.^3;
    Energy(N) = sum(P(:));
    P = Pold;
    RHO = 0.25*Pold.^2;
    Pold = double(Pold);
    %----------------------------------------------------------------------
    P = double(P);
    for k = 1:nb_sources
        tic;
            % P (nx,ny) --> Ob (nx+2,ny+2)
            % true/1 pour l'obstacle
            [U, Grad] = fm2dGradient(h, P, w, source_points(:,k), destination_points, Ob);
        toc
        for l=1:nb_dest
            dest = destination_points(:,l);
            dest = dest(1) + Nx*dest(2);
            Energy(N) = Energy(N) - Traffic(k,l)*U(dest);
            RHO = RHO - Traffic(k,l).*(squeeze(Grad(:,:,l)));
        end;
    end;
    P = P - rho*RHO;
    P = max(P,0);
    %----------------------------------------------------------------------
    if N > 1
        if Energy(N) > Energy(N-1)
            Nvar = Nvar + 1;
        end;
        if Nvar > 5
            K = K+1;
            rho = 1/K;
            Nvar = 0;
            ['At iteration ' int2str(N) ' \rho = ' num2str(rho) ';']
        end;
    end;
    %----------------------------------------------------------------------
    Pold = P;
    save(['MetricCavern/Ksi' int2str(N) '.mat'], 'P', 'Grad');
    imwrite(uint8(255*P/mmax(P)), ['MetricCavern/Ksi' int2str(N) '.png'], 'PNG');
    N = N+1;
    %----------------------------------------------------------------------
    figure(1);
    imshow(P,[min(P(:)) max(P(:))],'InitialMagnification',300);    
    hold on;
    plot(source_points(2,:),source_points(1,:),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',DisplaySourcesColor,'MarkerSize',DisplayMarkerSize);
    plot(destination_points(2,:),destination_points(1,:),'o',...
        'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',DisplayMarkerSize);
    colormap(jet);
    colorbar('vertical');
    axis on;
    set(gca,'XTickLabel',{''})
    set(gca,'XTick',[]);
    set(gca,'YTickLabel',{''})
    set(gca,'YTick',[]);
    title(['iter ',num2str(N-1) 'P']);
    pause(0.01);
    %----------------------------------------------------------------------
end;
%%
figure(3)
plot(Energy);
title('Energy');

%%
[X Y] = meshgrid(1:Nx, 1:Ny);
figure; contour3(X, Y, double(P), 40)
hold on
surf(X, Y, double(P), 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9])
axis off
set(gcf, 'Color', 'w')
colormap(gray(256));
light;

%% square logarithmic colormap
CC  = 1:(exp(1)-1)/1024:exp(1);
CM = log(CC); CM = [CM;CM;CM];
load 'MetricCavern/Ksi399.mat';
P(32,188) = mmax(P);
Obb = Ob(2:Nx+1, 2:Ny+1);
P(find(Obb)) =NaN;
figure;
set(gcf,'color',[1 1 1]);
imshow(log(P), [mmin(log(P)) mmax(log(P))],'InitialMagnification', 100);
colormap(1-sqrt(CM'));
hold on; contour(P, 31, 'Color', [0 0 0.01]);
colorbar;
hold on;
contour(double(Obb), 1,'color', [0 0 0.01]);
hold on;
axis on;
set(gca,'XTickLabel',{''})
set(gca,'XTick',[]);
set(gca,'YTickLabel',{''})
set(gca,'YTick',[]);
colorbar('vertical');
%% square logarithmic colormap
%CC  = 1:(exp(1)-1)/1024:exp(1);
%CM = log(CC); 
CM = rand(64,1);
CM = [CM';CM';CM'];
load 'MetricCavern/Ksi399.mat';
P(32,188) = mmax(P);
Obb = Ob(2:Nx+1, 2:Ny+1);
P(find(Obb)) =NaN;
figure;
set(gcf,'color',[1 1 1]);
imshow(log(P), [mmin(log(P)) mmax(log(P))],'InitialMagnification', 100);
colormap(1-CM');
hold on; contour(P, 31, 'Color', [0 0 0.01]);
colorbar;
hold on;
contour(double(Obb), 1,'color', [0 0 0.01]);
hold on;
axis on;
set(gca,'XTickLabel',{''})
set(gca,'XTick',[]);
set(gca,'YTickLabel',{''})
set(gca,'YTick',[]);
colorbar('vertical');


%%
for i=1:360
 camorbit(1,0,'data',[0 0 1])
 drawnow;
 %rotate(h, [0 0 1], 1);
 pause(0.01);
 hold on
 frame = getframe(gcf);
 imap = gray(256);   
 if i==1
        imwrite(squeeze(frame.cdata(:,:,1)),imap,  'RotateGray.gif','gif',...
            'DelayTime',1/20, 'LoopCount',Inf,...
            'WriteMode','overwrite')
    else
        imwrite(squeeze(frame.cdata(:,:,1)), imap, 'RotateGray.gif','gif','WriteMode','append',...
            'DelayTime',1/20)
    end
 pause(0.01);
end;