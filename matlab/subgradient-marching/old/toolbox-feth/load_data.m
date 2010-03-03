clear all;
close all;
clc;
stack=tiffread2('biodata01_crop2.tif');
%==============================================================
sz=[stack(1).height stack(1).width length(stack)];
img=zeros(sz);
for z=1:sz(3)
    img(:,:,z)=single(stack(z).data);
end;
clear stack;
min_img=min(img(:));
img=img-min_img;
max_img=max(img(:));
img=img./max_img;
%==============================================================
im = img(:, : , 18:140);
min_im = min(im(:));
im = im - min_im;
max_im = max(im(:));
im = im./max_im;
%==============================================================

%SliceIndex=50;
ù%figure;
%set(gcf,'color',ones(1,3));
%imshow(im(:,:,SliceIndex),[0 1],'InitialMagnification','fit');
%colormap(gray(256));
%axis on;
%set(gca,'XTickLabel',{''})
%set(gca,'XTick',[]);
%set(gca,'YTickLabel',{''})
%set(gca,'YTick',[]);