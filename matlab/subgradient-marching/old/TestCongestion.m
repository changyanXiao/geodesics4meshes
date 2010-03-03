%--------------------------------------------------------------------------
clear all; close all; clc;
%--------------------------------------------------------------------------
if 0
    library_directory = '../../fmlib/';
addpath(library_directory);
data_directory   = '../../Data/';
addpath(data_directory);
end
matlab_directory   = 'matlab/';
addpath(matlab_directory);
%--------------------------------------------------------------------------
DisplaySourcesColor=[1 0 0];
DisplayKeyPointsColor=[1 1 0];
DisplayMarkerSize=4;
%--------------------------------------------------------------------------
Nx = 101; Ny = 101;
Ob = zeros(Nx+2, Ny+2);
Ob = logical(Ob);
Pold = ones(Nx,Ny);
%--------------------------------------------------------------------------
% animated gif options
%FileName = 'CongAnim.gif';
ColorMap = jet(256);
framesPerSecond = 20;
loopCount = Inf;
%--------------------------------------------------------------------------
w = 0.001;
h = [1;1];
source_points = [51;15];
destination_points = [51;87];
sz_source = size(source_points); nb_sources = sz_source(2);
sz_dest   = size(destination_points); nb_dest = sz_dest(2);
Traffic = ones(nb_sources, nb_dest);
Nb_iter = 1000;
Energy = ones(1,Nb_iter-1);
N = 1;
Nvar = 0;
K = 5;
rho = 1/K;
%--------------------------------------------------------------------------
while(N < Nb_iter)
    %----------------------------------------------------------------------
%    if N==1
%        imwrite(uint8(255*Pold/mmax(Pold)),ColorMap,FileName, 'gif',...
%            'DelayTime',1/framesPerSecond, 'LoopCount',loopCount,...
%            'WriteMode','overwrite')
%    else
%        imwrite(uint8(255*Pold/mmax(Pold)),ColorMap,FileName,'gif','WriteMode','append',...
%            'DelayTime',1/framesPerSecond)
%    end;
    %----------------------------------------------------------------------
    P = (1.0/3.0).*Pold.^3;
    Energy(N) = sum(P(:));
    P = Pold;
    RHO = Pold.^2;
    Pold = double(Pold);
    %----------------------------------------------------------------------
    for k = 1:nb_sources
        for l = 1:nb_dest
            P = double(P);
            tic;
                [U, Grad] = fm2dGradient(h, P, w, source_points, destination_points, Ob);
            toc
            Grad = squeeze(Grad);
            Energy(N) = Energy(N) - Traffic(k,l)*U(51,87);
            RHO = RHO - Traffic(k,l).*(Grad);
            P = P - rho*RHO;
            P = max(P,0);
        end;
    end;
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
    save(['Metric2/Ksi' int2str(N) '.mat'], 'P', 'Grad');
    imwrite(uint8(255*P/mmax(P)), ['Metric2/Ksi' int2str(N) '.png'], 'PNG');
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
plot(Energy(1:300));
title('Energy');
%% plot 3D elevated surface
P(51,15) = mmax(P);
[X Y] = meshgrid(1:Nx, 1:Ny);
figure; contour3(X, Y, double(P), 30)
hold on
surf(X, Y, double(P), 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9])
axis off
set(gcf, 'Color', 'w')
colormap(gray(256));
light;
%% Rotate and save animated gif
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
%% square logarithmic colormap
CC  = 1:(exp(1)-1)/256:exp(1);
CM = log(CC); CM = [CM;CM;CM];
figure;
set(gcf,'color',[1 1 1]);
imshow(P, [mmin(P) mmax(P)]);
colormap(1-sqrt(CM'));
hold on; contour(P, 30, 'Color', [0 0 0.01]);
colorbar;
%% random colormap
CC  = rand(1,64);
%CM =(1./CC);
CM = [CC;CC;CC];
figure; imshow(P, [mmin(P) mmax(P)]);
colormap(1-CM')
hold on; contour(P, 30, 'Color', [0 0 0.01]);
colorbar;
set(gcf, 'Color', 'w');

%% logarithmic colormap with thresholded matric
CC  = 1:(exp(1)-1)/255:exp(1);
CM = log(CC); CM = [CM;CM;CM];
PP = P;
PP(P > mmax(P)/3) = mmax(P)/3;
figure; imshow(PP, [mmin(PP) mmax(PP)]);
colormap(1-CM')
hold on; contour(PP, 30, 'Color', [0 0 0.01]);
colorbar;
set(gcf, 'Color', 'w');
