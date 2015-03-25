%thresholding of peak duration and peak height

% clear all;clc;close all;
% run setup_header.m;
% rootpath='C:\nuclei\post analysis result_0.2\data';
% fileswt=dir([rootpath,'\wild_type*.mat']);
% for ifile=1%:length(fileswt)
%     load(fullfile(rootpath,fileswt(ifile).name));
%
% end
clc
clear all
run setup_header.m;
rootpath='C:\nuclei\post analysis result_0.2';
load([rootpath,'\strainall.mat']);

%%
phth=0.05;
durth=25;
clf
fluc=strainall(1).flucs;
fluc=fluc([fluc.good]==1);
initt=zeros(size(fluc));
for i=1:length(fluc)
    initt(i)=sum(fluc(i).width~=0);
end
plot([fluc.meanheight],initt*2.5,'r.','markersize',20);hold on;
%     plot([fluc.maxheight],[fluc.risetime]+[fluc.falltime],'g.');hold on;
wt_out_ratio=sum(initt*2.5>=25 &[fluc.meanheight]>0.05)/length(initt)



fluc=strainall(9).flucs;
fluc=fluc([fluc.good]==1);
initt=zeros(size(fluc));
for i=1:length(fluc)
    initt(i)=sum(fluc(i).width~=0);
end
plot([fluc.meanheight],initt*2.5,'.','markersize',20,'color',[1 1 1]*0);hold on;
%     plot([fluc.maxheight],[fluc.risetime]+[fluc.falltime],'b.');
sp10mbc_out_ratio=sum(initt*2.5>=25&[fluc.meanheight]>0.05)/length(initt)

plot([.05 .05], [25 100],'b');hold on;
plot([.05 .2], [25 25],'b');hold on;

axis([0 .15 0 100])
legend('wt','sp10 MBC');
xlabel('flucutation height ( \mum )');
ylabel('duration of fluctuation above threshold ( s )');
set(gca,'YTick',0:20:100)
FigureFormat(gcf)