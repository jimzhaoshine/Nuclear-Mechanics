% plot pictures for sarah
% wild type certain frames

load('C:\nuclei\post analysis result_0.2\data\wild_type_04.mat');
% nm.path=['C',nm.path(2:9),'\data',nm.path(10:end)];
nm.loadmovie;
%%
inuc=3;
iframe=18;
%% 3d
for iframe=1:1:30;
    clf
    nuc=nm.nuclei{iframe,inuc};
    cnt=round(nuc.origin);
    points=nm.points;
    faces=nm.faces;
    pts=[nuc.r.*points(:,1),-nuc.r.*points(:,2),nuc.r.*points(:,3)];
    TR = triangulation(faces,pts);
    trisurf(TR,'FaceColor','red','EdgeColor','none');
    axis([-15 15 -15 15 -15 15]);
%     axis off
    % view([-0.5 -1 0]);
    view([0 0 1]);
    daspect([1 1 1/0.85])
    camlight
    lighting gouraud
    savefig(['pic1\3dwtflucframe',num2str(iframe)]);
end
%                                 set(ahandle,'box','on','XColor',[1 1 1],'YColor',[1 1 1])
%                                 set(ahandle,'xtick',[],'ytick',[],'linewidth',1)

%% 2d
wsize=25;
for iframe=1:1:30;
    img=nm.grab3(iframe);
    nuc=nm.nuclei{iframe,inuc};
    cnt=round(nuc.origin);
    %     istack=cnt(3)-1:cnt(3)+1;
    istack=cnt(3)-1:cnt(3)+1;
    for i=1:3
    dx=nuc.origin(1)-round(nuc.origin(1));
    dy=nuc.origin(2)-round(nuc.origin(2));
    if istack(i)<=10 && istack(i)>0
        pos=[nuc.contour(istack(i)).x+dx,nuc.contour(istack(i)).y+dy]+(wsize+1);
        %         ahandle=axes('Units','Pixels','Position',[(isub_nuc-1)*sfw i*sfw-sfw sfw sfw]);
        wimg=WindowImageUS(cnt(1),cnt(2),wsize,squeeze(img(:,:,istack(i))));
        imagesc(wimg);axis image;colormap gray;hold on;
        plot(pos(:,1),pos(:,2),'r-','linewidth',2);    axis off;
        
        surf(wimg'/max(wimg(:)));view([1 0 5]);axis off;
        
        
    end
    print(['pic1\wtfluc_frame',num2str(iframe),'_stack',num2str(istack(i))],'-dtiff');
%     pause
    end
end

%% kemograph
extsiz=10;
% half window size for particle finding, centroid caculation,
% periodic expansion
hpsiz=8;
% calculate mean contour
sum_r=0;
for iframe=1:nm.endframe
    sum_r=sum_r+nm.nuclei{iframe,inuc}.r_new;
end
mean_r=sum_r/nm.endframe;
meannuc=nm.nuclei{1,inuc};
meannuc.r_new=mean_r;
[ points2,faces2,psi,theta ] = setup_img(  );
meannuc=centralband(meannuc,nm,points2,size(psi),1);
mean_img=meannuc.img;
% construct 3d image series and filtered image;
time_img=zeros(size(psi,1)+2*hpsiz,size(psi,2)+2*extsiz+2*hpsiz,nm.endframe);
for iframe=1:nm.endframe
    nuc=nm.nuclei{iframe,inuc};
    img_filtered=nuc.img;
    dimg=(img_filtered-mean_img)./mean_img;
    dimg=[dimg(:,end-hpsiz+1:end),dimg,dimg(:,1:hpsiz)];
    dimg=[zeros(hpsiz,size(dimg,2));...
        dimg;zeros(hpsiz,size(dimg,2))];
    time_img(:,:,iframe)=[dimg(:,end-extsiz+1:end),dimg,dimg(:,1:extsiz)];
end
% filter the data
time_img_filtered=gausspass3(time_img,2,2);
time_img_filtered=time_img_filtered(:,extsiz+1:end-extsiz,:);
time_img=time_img(:,extsiz+1:end-extsiz,:);
% calculate std
std_img=std(time_img,1,3);
%%
figure('Position',[0 0 600 800]);
img_centerline=(squeeze(time_img_filtered(10,:,:)))';
%             img_centerline=[img_centerline(:,end-hpsiz+1:end),img_centerline,img_centerline(:,1:hpsiz)];
imagesc(img_centerline,[-0.05,0.2]);colormap jet;
xlabel('longitudinal angle')
ylabel('time (s)')
set(gca,'XTick',(0:16:64)+hpsiz);
set(gca,'YTick',(0:20:100),'YTicklabel',(0:20:100)*2.5);
title('longitutde centerline fluctuation kymograph')
set(gca, 'xticklabel', '0 | p/2 | p | 3p/2 | 2p', 'fontname', 'symbol');
colorbar;
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(gcf,'type','axes'),'fontSize',14,'fontWeight','bold')
print('pic1\centerkemo','-deps');
print('pic1\centerkemo','-dpng');
savefig('pic1\centerkemo');
%%
figure('Position',[0 0 600 800]);
max_longi=(squeeze(max(time_img_filtered,[],1)))';
%             max_longi=[max_longi(:,end-hpsiz+1:end),max_longi,max_longi(:,1:hpsiz)];
imagesc(max_longi,[-0.05,0.2]);colormap jet;hold on;
xlabel('longitudinal angle')
ylabel('time (s)')
set(gca,'XTick',(0:16:64)+hpsiz);
set(gca,'YTick',(0:20:100),'YTicklabel',(0:20:100)*2.5);
title('longitutde max fluctuation kymograph')
set(gca, 'xticklabel', '0 | p/2 | p | 3p/2 | 2p', 'fontname', 'symbol');
colorbar
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(gcf,'type','axes'),'fontSize',14,'fontWeight','bold')
print('pic1\longitudemaxkemo','-deps');
print('pic1\longitudemaxkemo','-dpng');
savefig('pic1\longitudemaxkemo');
%%
figure('Position',[0 0 1200 600]);
std_img1=std_img(1+hpsiz:end-hpsiz,:);
imagesc(std_img1,[0.02 0.1]);colormap jet;axis image;colorbar;
title('standard deviation map');hold on;
xlabel('longitudinal angle')
ylabel('latitudinal angle')
set(gca,'XTick',(0:16:64));
set(gca, 'xticklabel', '0 | p/2 | p | 3p/2 | 2p', 'fontname', 'symbol');
set(gca,'YTick',([2,10,18]));
set(gca, 'yticklabel', '-p/4 | 0 | p/4', 'fontname', 'symbol');
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(gcf,'type','axes'),'fontSize',14,'fontWeight','bold')

print('pic1\stdimg','-deps');
print('pic1\stdimg','-dpng');
savefig('pic1\stdimg');
