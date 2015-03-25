%% time coorelation analysis
z_range=5;
points=nm.points;
faces=nm.faces;
neighbors=nm.neighbors;
[ points2,faces2,psi,theta] = setup_img( );
pix2um=0.16;
% calculate mean 3d contour
for inuc=1:nm.num_nuc
    %% calculate mean contour
    sum_r=0;
    for iframe=1:nm.endframe
        sum_r=sum_r+nm.nuclei{iframe,inuc}.r_new;
    end
    mean_r=sum_r/nm.endframe;
    mean_pts=[mean_r.*points(:,1),mean_r.*points(:,2),mean_r.*points(:,3)];
    patchm.vertices=mean_pts;
    patchm.faces=faces;
    meannuc=nm.nuclei{1,inuc};
    meannuc.r_new=mean_r;
    meannuc=centralband(meannuc,nm,points2,size(psi));
    mean_img=meannuc.img;
    
    %  check center position
    for vt=[]
        center_smooth=[];
        center_noisy=[];
        for iframe=1:nm.endframe
            center_smooth(iframe,:)=nm.nuclei{iframe,inuc}.origin_new;
            center_noisy(iframe,:)=nm.nuclei{iframe,inuc}.center;
        end
        C=cov(center_smooth(:,1),center_smooth(:,2));
        [V,D]=eig(C);
        [~,vind]=max([D(1),D(4)]);
        eigv=V(:,vind);
        mx=mean(center_smooth(:,1));
        my=mean(center_smooth(:,2));
        mz=mean(center_smooth(:,3));
        nuc_psi=acos(eigv(1));
        if 0
            plot3(center_smooth(:,1),center_smooth(:,2),center_smooth(:,3),'bo-');hold on;
            plot3(center_noisy(:,1),center_noisy(:,2),center_noisy(:,3),'ro-');hold on;
            plot3(mx+[-1,1]*eigv(1),my+[-1,1]*eigv(2),[0 0]+mz,'--g')
            legend('filtered','raw centroid')
            title('drift control filtering');
        end
        if 0
            figure(2)
            clf
            plot(center_smooth(:,1),center_smooth(:,2),'bo-');hold on;
            plot(center_noisy(:,1),center_noisy(:,2),'ro-');hold on;
            plot(mx+[-1,1]*eigv(1),my+[-1,1]*eigv(2),'--g');
            legend('filtered','raw centroid')
            title('drift control filtering');
        end
    end
    
    % 2d analysis
    time_img=zeros(size(psi,1),size(psi,2),nm.endframe);
    abs_time_img=zeros(size(psi,1),size(psi,2),nm.endframe);
    for iframe=1:nm.endframe
        nuc=nm.nuclei{iframe,inuc};
        img=nuc.img;
        dimg=img-mean_img;
        time_img(:,:,iframe)=dimg./mean_img;
        abs_time_img(:,:,iframe)=img;
    end
    % calculate std
    std_img=std(time_img,1,3);
    %plot std
    for vt=[]
    f=figure(1);
    clf
    set(f,'Position',[100 100 1200 600]);
    %     axes('Units','Pixels','Position',[50 50 ]);
    subplot(2,1,1)
    hpsiz=4;
    ext_std_img=[std_img(:,end-hpsiz+1:end),std_img,std_img(:,1:hpsiz)];
    ext_std_img=[zeros(hpsiz,size(time_img,2)+2*hpsiz)+min(std_img(:));...
        ext_std_img;zeros(hpsiz,size(time_img,2)+2*hpsiz)+min(std_img(:))];
    
    %peak finding
    pks=pkfnd(ext_std_img,0.04,hpsiz);
    pks=pks-hpsiz;
    
    std_img=std_img(5:end-4,:);
    imagesc(std_img,[0.02, .1]);colormap jet;axis image;hold on;
    plot(nuc_psi/2/pi*64+[0 0],[0,size(std_img,1)+1],'r-');hold on;
    plot(nuc_psi/2/pi*64+32+[0 0],[0,size(std_img,1)+1],'r-');hold on;
    if ~isempty(pks)
        pks(:,2)=pks(:,2)-4;
        pks=pks(pks(:,1)>0&pks(:,2)>0&pks(:,1)<=size(std_img,2)&pks(:,2)<=size(std_img,1),:);
        plot(pks(:,1),pks(:,2),'ko');
    end
    xlabel('\psi');
    ylabel('\theta');
    set(gca,'XTick',0:16:65);
    set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
    set(gca,'YTick',[1 6 11]);
    set(gca,'YTickLabel',{'5pi/32','0','-5pi/32'})
    colorbar;
    title('standard deviation');
    %     axes('Position',[0 0 1 1/2]);
    
    subplot(2,1,2)
    if ~isempty(pks)
        flucs=[];
        for i=1:size(pks,1)
            flucs(:,i)=squeeze(squeeze(time_img(pks(i,2),pks(i,1),:)));
        end
        plot((1:nm.endframe)*2.5,flucs);hold on;
        plot([0 nm.endframe*2.5],[0 0 ],'--g')
        xlabel('times (s)');
        ylabel('relative fluctuation (dr/mean)');
    end
    end

    %% filter data and plot kemo
img3=gausspass3(time_img,1,3);
clf;subplot(1,2,1)
imagesc(squeeze(img3(10,:,:))',[-0.05,0.2]);colormap jet;
xlabel('angle');ylabel('frames'); set(gca,'XTick',10:16:74);
set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
title('longitutde centerline fluctuation kemograph');colorbar
subplot(1,2,2);max_longi=squeeze(max(img3,[],1));
imagesc(max_longi',[-0.05,0.2]);colormap jet;
xlabel('angle');ylabel('frames');set(gca,'XTick',10:16:74);
set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
title('longitutde max fluctuation kemograph');colorbar

    %% track peaks for each image
    pos=[];
    for iframe=1:nm.endframe
        img=squeeze(img3(:,:,iframe));
        hpsiz=4;
        ext_img=[img(:,end-hpsiz+1:end),img,img(:,1:hpsiz)];
        ext_img=[zeros(hpsiz,size(ext_img,2))+mean(img(:));...
            ext_img;zeros(hpsiz,size(ext_img,2))+mean(img(:))];
        pks=pkfnd(ext_img,0.07,hpsiz);
        cnt=cntrd(ext_img,pks,hpsiz);
        pos=[pos;cnt,iframe+zeros(size(cnt,1),1)];
    end
    maxdisp=3;
    tparam.mem=2;
    tparam.dim=2;
    tparam.good=5;
    tparam.quiet=0;
    tracks=track(pos(:,[1 2 5]),maxdisp,tparam);
    
    max_longi=squeeze(max(img3,[],1));
    imagesc(max_longi',[-0.05,0.2]);colormap gray;hold on;
    xlabel('angle');ylabel('frames');set(gca,'XTick',10:16:74);
    set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
    title('longitutde max fluctuation kemograph');colorbar
     for itj=1:tracks(end,4)
        trtmp=tracks(tracks(:,4)==itj,:);
        xtr=trtmp(:,1)-hpsiz;
        ytr=trtmp(:,2)-hpsiz;
        ttr=trtmp(:,3);
        plot(xtr,ttr,'o-','Color',GenColor(itj/tracks(end,4)));hold on;
    end
    
   %     SI(std_img);hold on;
%  for itj=1:tracks(end,4)
%         trtmp=tracks(tracks(:,4)==itj,:);
%         xtr=trtmp(:,1)-hpsiz;
%         ytr=trtmp(:,2)-hpsiz;
%         plot(xtr,ytr,'o-','Color',GenColor(itj/tracks(end,4)));hold on;
%     end
    
    %% plot 2d fluctuation over time
    for vt=[]
    f=figure(3);
    set(f,'Position',[100 100 1300 500]);
    for iframe=1:nm.endframe
        nuc=nm.nuclei{iframe,inuc};
        img=nuc.img;
        relative_fluc_img=(img-mean_img)./mean_img;
        relative_fluc_img=relative_fluc_img(5:end-4,:);
        clf
        a1=axes('Position',[0 0.05 1 0.4]);
        imagesc(relative_fluc_img,[-0.2, .2]);colormap jet;axis image;hold on;
        colorbar;
        set(gca,'XTick',0:16:65);
        set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
        set(gca,'YTick',1:5:11);
        set(gca,'YTickLabel',{'5pi/32','0','-5pi/32'})
        title(['relative fluctuation at frame ',num2str(iframe)])
        
        a2=axes('Position',[0 0.55 1 0.4]);
        imagesc(std_img,[0.0, .1]);colormap jet;axis image;hold on;
        colorbar;
        set(gca,'XTick',0:16:65);
        set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
        set(gca,'YTick',[1 6 11]);
        set(gca,'YTickLabel',{'5pi/32','0','-5pi/32'})
        title('relative fluctuation standard deviation')
        pause
    end
    end
end

%% make movie for each 2d band
daObj=VideoWriter(['centralband ',nm.filename,' ',num2str(inuc),'th'],'MPEG-4');
daObj.set('Quality',100);
daObj.FrameRate=1;
open(daObj);

f=figure(3);
set(f,'Position',[100 100 1300 500]);
for iframe=1:nm.endframe
    
    nuc=nm.nuclei{iframe,inuc};
    img=nuc.img;
    relative_fluc_img=(img-mean_img)./mean_img;
    relative_fluc_img=relative_fluc_img(5:end-4,:);
    clf
    a1=axes('Position',[0 0.05 1 0.4]);
    imagesc(relative_fluc_img,[-0.2, .2]);colormap jet;axis image;hold on;
    colorbar;
    %         xlabel('\psi');
    %     ylabel('\theta');
    set(gca,'XTick',0:16:65);
    set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
    set(gca,'YTick',1:5:11);
    set(gca,'YTickLabel',{'5pi/32','0','-5pi/32'})
    title(['relative fluctuation at frame ',num2str(iframe)])
    
    a2=axes('Position',[0 0.55 1 0.4]);
    imagesc(std_img,[0.0, .1]);colormap jet;axis image;hold on;
    colorbar;
    %         xlabel('\psi');
    %     ylabel('\theta');
    set(gca,'XTick',0:16:65);
    set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
    set(gca,'YTick',[1 6 11]);
    set(gca,'YTickLabel',{'5pi/32','0','-5pi/32'})
    title('relative fluctuation standard deviation')
    writeVideo(daObj,getframe(f));
end
close all;
close(daObj);



%%

% max fluc at diff longitude at time series
max_longi=squeeze(max(time_img,[],1));
SI(max_longi');colormap jet;axis image
xlabel('angle')
ylabel('frames')
set(gca,'XTick',10:16:74);
set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
title('longitutde max fluctuation time series')
colorbar

%% 2d particle tracking
time_img=zeros(17,84,101);
for iframe=1:nm.endframe
    nuc=nm.nuclei{iframe,1};
    img=nuc.img;
    dimg=img-mean_img;
    dimg=[dimg(:,end-1:end),dimg,dimg(:,1:2)];
    pks=pkfnd(dimg,.2,2);
    imagesc(dimg,[-1 1]);colormap jet; axis image;
    hold on;
    if ~isempty(pks)
        plot(pks(:,1),pks(:,2),'ko')
    end
    pause;
end
%% pks group data
imgg=[];
for i=1:size(pks,1)
    img1=squeeze(time_img(pks(i,2),pks(i,1)+ (-10:10),:))';
    imgg=[imgg,img1];
end
SI(imgg);colormap jet;
%% 3d analysis
for inuc=1:nm.num_nuc
    f=figure('Position',[100 100 1200 600]);
    for iframe=1:nm.endframe
        
        nuc=nm.nuclei{iframe,inuc};
        clf
        dr=nuc.r_new-mean_r;
        r=nuc.r_new;
        x=r.*points(:,1);
        y=r.*points(:,2);
        z=r.*points(:,3);
        cntall=[x,y,z];
        
        % find pks of fluctuations in 3d
        pks=[];
        pksi=[];
        ind_center=find(z>-z_range & z<z_range);
        for indi=1:length(ind_center)%1:size(r)
            i=ind_center(indi);
            ri=r(i);
            dri=dr(i);
            dri1max=max(dr(neighbors(i,~isnan(neighbors(i,:)))));
            %             ri=r(i);
            %             cnt=cntall(i,:);
            %             ri1=mean(r(neighbors(i,~isnan(neighbors(i,:)))));
            %             cnt1=mean(cntall(neighbors(i,~isnan(neighbors(i,:))),:));
            %             ri1m=max(r(neighbors(i,~isnan(neighbors(i,:)))));
            %             ri2=mean(r(neighbors2(i,~isnan(neighbors2(i,:)))));
            %             cnt2=mean(cntall(neighbors2(i,~isnan(neighbors2(i,:))),:));
            %             d12=sqrt((cnt-cnt1)*(cnt1-cnt2)');
            if  dri>0.2 && dri>=dri1max
                pks=[pks;ri*points(i,1),ri*points(i,2),ri*points(i,3)];
                %                 plot3(ri*points(i,1),ri*points(i,2),ri*points(i,3),'ob');hold on;
            end
        end
        
        % plot 3d nuclei and mean contour
        figure(f)
        a1=axes('Units','Pixels','Position',[0 0 600 600]);
        pts1=[r.*points(:,1),r.*points(:,2),r.*points(:,3)];
        patch1.vertices=pts1;
        patch1.faces=faces;
        pts2=[mean_r.*points(:,1),mean_r.*points(:,2),mean_r.*points(:,3)];
        %         pts2=pts2+ones(size(pts2,1),1)*corcnt(iframe,:);
        patch2.vertices=pts2;
        patch2.faces=faces;
        %
        patch(patch1,'FaceColor','red','EdgeColor','none','FaceAlpha',.3);hold on;
        patch(patch2,'FaceColor','blue','EdgeColor','none','FaceAlpha',.3);hold on;
        if ~isempty(pks)
            plot3(pks(:,1),pks(:,2),pks(:,3),'og');hold on;
        end
        plot3([-10 10],[0 0],[z_range z_range],'k');hold on;
        plot3([-10 10],[0 0],[-z_range -z_range],'k');hold on;
        
        view([1 1 0]);
        axis([-10 10 -10 10 -10 10])
        daspect([1 1 1])
        grid off
        axis off
        camlight
        lighting gouraud
        
        a2=axes('Units','Pixels','Position',[600 0 600 600]);
        pts1=[r.*points(:,1),r.*points(:,2),r.*points(:,3)];
        patch1.vertices=pts1;
        patch1.faces=faces;
        pts2=[mean_r.*points(:,1),mean_r.*points(:,2),mean_r.*points(:,3)];
        %         pts2=pts2+ones(size(pts2,1),1)*corcnt(iframe,:);
        patch2.vertices=pts2;
        patch2.faces=faces;
        %
        patch(patch1,'FaceColor','red','EdgeColor','none','FaceAlpha',.3);hold on;
        patch(patch2,'FaceColor','blue','EdgeColor','none','FaceAlpha',.3);hold on;
        if ~isempty(pks)
            plot3(pks(:,1),pks(:,2),pks(:,3),'og');hold on;
        end
        
        view([0 0 1]);
        camzoom(1/sqrt(2))
        axis([-10 10 -10 10 -10 10])
        daspect([1 1 1])
        grid off
        axis off
        camlight
        lighting gouraud
        
        pause
    end
end
