% plot pictures for sarah
% make all pictures and save
clear all;
close all;
clc;
savedir='C:\nuclei\post analysis result_0.2\save';
mkdir(savedir);
%%
% path='C:\nuclei\post analysis result_0.2\data\cdc25';
path='C:\nuclei\post analysis result_0.2\data\cdc25\';
files=dir([path,'*.mat']);
for ifile=1:length(files)
moviename=files(ifile).name(1:end-4);
load([path,moviename,'.mat']);
nm.loadmovie(10);
iframe=1;
for inuc=1:nm.num_nuc
%% 3d
clf
nuc=nm.nuclei{iframe,inuc};
cnt=round(nuc.origin);
points=nm.points;
faces=nm.faces;
pts=[nuc.r.*points(:,1),-nuc.r.*points(:,2),nuc.r.*points(:,3)];
TR = triangulation(faces,pts);
trisurf(TR,'FaceColor','red','EdgeColor','none');
axis([-15 15 -15 15 -15 15]);
view([0 0 1]);
daspect([1 1 1/0.85])
% axis off;
grid off;
set(gca,'Xtick',[-12.5 -6.25 0 6.25 12.5],'Ytick',[-12.5 -6.25 0 6.25 12.5]);
set(gca, 'xticklabel', [{'-2'}, {-1},{'0'},{'1'},{'2'}], 'yticklabel', [{'-2'}, {-1},{'0'},{'1'},{'2'}] );
xlabel('\mum')
ylabel('\mum')
camlight
lighting gouraud
% print([savedir,'\3d_',moviename,'_',num2str(inuc),'.png'],'-dpng');
savefig([savedir,'\3d_',moviename,'_',num2str(inuc),'.fig']);

%% 2d
wsize=20;
img=nm.grab3(iframe);
nuc=nm.nuclei{iframe,inuc};
cnt=round(nuc.origin);
dx=nuc.origin(1)-round(nuc.origin(1));
dy=nuc.origin(2)-round(nuc.origin(2));
for ii=1:3
    
istack=round(cnt(3))+(ii-2);
if istack>=1 && istack<=10
pos=[nuc.contour(istack).x+dx,nuc.contour(istack).y+dy]+(wsize+1);
wimg=WindowImageUS(cnt(1),cnt(2),wsize,squeeze(img(:,:,istack)));
imagesc(wimg);axis image;colormap gray;hold on;
plot(pos(:,1),pos(:,2),'r-','linewidth',2);   
set(gca,'Xtick',[-12.5 -6.25 0 6.25 12.5]+wsize+1,'Ytick',[-12.5 -6.25 0 6.25 12.5]+wsize+1);
set(gca, 'xticklabel', [{'-2'}, {-1},{'0'},{'1'},{'2'}], 'yticklabel',[{'-2'}, {-1},{'0'},{'1'},{'2'}] );
xlabel('\mum')
ylabel('\mum')

% print([savedir,'\2d_',moviename,'_',num2str(inuc),'.png'],'-dpng');
savefig([savedir,'\2d_',moviename,'_',num2str(inuc),'_',num2str(ii),'.fig']);
end
end
end
end