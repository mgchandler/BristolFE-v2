addpath([pwd '\mtex-5.2.3' ]);startup;
% Randomly generate a set of grain orientations and plot pole figure and inverse pole figure 
cs = crystalSymmetry('cubic');
ss=specimenSymmetry('triclinic');
N=100;
ori = orientation.rand(N,cs,ss);
figure;
plotPDF(ori,Miller({1,0,0},{1,1,0},{1,1,1},cs),'all','smooth');
% plotPDF(ori1,Miller({1,0,0},{1,1,0},{1,1,1},cs),'all');
mtexColorbar;
vec=[vector3d.X,vector3d.Y,vector3d.Z];
figure;
% plotIPDF(ori1,vec,'all','smooth');
plotIPDF(ori,vec,'all');
mtexColorbar;


