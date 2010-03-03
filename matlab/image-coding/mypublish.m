% publish files

opts.outputDir = 'html/';
opts.format = 'latex';
file = publish('test_geometric_coder',opts);

cd('html');
!latex test_geometric_coder
!dvips test_geometric_coder.dvi -o  test_geometric_coder.ps
!ps2pdf test_geometric_coder.ps
cd('..');