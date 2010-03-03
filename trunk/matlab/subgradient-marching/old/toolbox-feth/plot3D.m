function plot3D(Image, alpha)

h = vol3d('cdata', Image, 'texture', '3D');
view(3); 
vol3d(h);
axis tight;  %daspect([1 1 1]);
alphamap('rampup');
alphamap(alpha .* alphamap);
