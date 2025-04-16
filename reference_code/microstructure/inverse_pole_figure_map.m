addpath([pwd '\mtex-5.2.3' ]);startup;
% Randomly assign grain orientation to each grain of a known polycrystal ,and plot inverse pole figure map
load('grainid.mat')
cs=crystalSymmetry('cubic');
ss=specimenSymmetry('triclinic');
ipfkey=ipfColorKey(cs);
figure;
plot(ipfkey);
num_id=max(unique(grainid));
ori=orientation.rand(num_id,cs,ss);
color=ipfkey.orientation2color(ori);
ipfmap_1=grainid;ipfmap_2=grainid;ipfmap_3=grainid;
for i=1: num_id    
    ipfmap_1(ipfmap_1==i)=color(i,1);
    ipfmap_2(ipfmap_2==i)=color(i,2);
    ipfmap_3(ipfmap_3==i)=color(i,3);
end
ipfmap(:,:,1)=ipfmap_1;
ipfmap(:,:,2)=ipfmap_2;
ipfmap(:,:,3)=ipfmap_3;
figure;imagesc(ipfmap);axis equal;