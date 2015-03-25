% make all the kemographs for max longtitude and equator line, save to
% fig

path='C:\nuclei\post analysis result_0.05';
files=dir(path);
f=figure('Position',[100 100 1200 900]);
nucind=cell(length(files)-2,1);
trj=cell(length(files)-2,1);
for ifile=3:length(files);
    filename=files(ifile).name;
    [~,name,ext]=fileparts(filename);
    display(name);
    if strcmp(ext,'.mat')
        load(fullfile(path,filename));
        if ifile==3
            [ points2,faces2,psi,theta] = setup_img( );
        end
       
        for inuc=1:nm.num_nuc
            %% calculate mean contour
            extsiz=3;
            hpsiz=8;
                            
            sum_r=0;
            for iframe=1:nm.endframe
                sum_r=sum_r+nm.nuclei{iframe,inuc}.r_new;
            end
            mean_r=sum_r/nm.endframe;
            meannuc=nm.nuclei{1,inuc};
            meannuc.r_new=mean_r;
            meannuc=centralband(meannuc,nm,points2,size(psi));
            mean_img=meannuc.img;
            time_img=zeros(size(psi,1),size(psi,2)+2*extsiz,nm.endframe);
            for iframe=1:nm.endframe
                nuc=nm.nuclei{iframe,inuc};
                img=nuc.img;
                dimg=(img-mean_img)./mean_img;
                time_img(:,:,iframe)=[dimg(:,end-extsiz+1:end),dimg,dimg(:,1:extsiz)];
            end
            % calculate std
            std_img=std(time_img,1,3);
            % filter the data
            img3=gausspass3(time_img,2,5);
            img3=img3(:,extsiz+1:end-extsiz,:);
%             time_img=time_img(:,extsiz+1:end-extsiz,:);
            % track peaks for each image
            pos=[];
            for iframe=1:nm.endframe
                img=squeeze(img3(:,:,iframe));
                ext_img=[img(:,end-hpsiz+1:end),img,img(:,1:hpsiz)];
                ext_img=[zeros(hpsiz,size(ext_img,2))+mean(img(:));...
                    ext_img;zeros(hpsiz,size(ext_img,2))+mean(img(:))];
%                 pks=pkfnd(ext_img,0.07,hpsiz);
                pks=pkfnd(ext_img,0.05,hpsiz);
                cnt=[];
                if ~isempty(pks)
                    for ipeak=1:size(pks,1)
                        wimg=ext_img(pks(ipeak,2)+(-hpsiz:hpsiz),pks(ipeak,1)+(-hpsiz:hpsiz));
                        bw=wimg(wimg>0.5*wimg(hpsiz+1,hpsiz+1));
                        [imy,imx]=meshgrid(pks(ipeak,2)+(-hpsiz:hpsiz),pks(ipeak,1)+(-hpsiz:hpsiz));
                        my=mean(mean(wimg.*imy))/mean(wimg(:));
                        mx=mean(mean(wimg.*imx))/mean(wimg(:));
                        cnt=[cnt;mx,my,wimg(hpsiz+1,hpsiz+1),sqrt(sum(bw(:)))];
                    end
%                     cnt=cntrd(ext_img,pks,hpsiz);
                    pos=[pos;cnt,iframe+zeros(size(cnt,1),1)];
                end
            end
%             maxdisp=4;
%             tparam.mem=0;
%             tparam.dim=2;
%             tparam.good=3;
%             tparam.quiet=0;
            tparam.dim=2;
            tparam.maxdisp=5;
            tparam.mem=2;
            tparam.good=10;

            tracks=[];
            tracksnew=[];
            if ~isempty(pos)
%                 try
%                     pos_new=[pos(:,1:2),pos(:,5),pos(:,3:4)];
%                     tracks=track(pos_new(:,1:3),maxdisp);%,tparam);
%                     [ tracksnew ] = MergeProperties(tracks,pos_new,2 );
                    ytracks=yaotrack(pos,tparam);
                    if ~isempty(ytracks)
                        tracks=ytracks(:,[1 2 end-1 end]);
                        tracksnew=ytracks;
                    end
%                 catch
%                 end
            end
            nucind{ifile-2}=[name,'_',num2str(inuc)];
            trj{ifile-2}=tracksnew;
            
            % plot
            clf
            axes('Position',[0.05 0.4 0.4 0.55]);
            img1=(squeeze(img3(10,:,:)))';
            img1=[img1(:,end-hpsiz+1:end),img1,img1(:,1:hpsiz)];
            imagesc(img1,[-0.05,0.2]);colormap jet;
            xlabel('angle')
            ylabel('frames')
            set(gca,'XTick',(1:16:64)+hpsiz);
            set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
            title('longitutde centerline fluctuation kemograph')
            colorbar;hold on;
            if ~isempty(tracks)
                for itj=1:tracks(end,4)
                    trtmp=tracks(tracks(:,4)==itj,:);
                    xtr=trtmp(:,1);%-hpsiz;
                    ytr=trtmp(:,2)-hpsiz;
                    ttr=trtmp(:,3);
                    plot(xtr,ttr,'.-k');hold on;%,'Color',GenColor(itj/tracks(end,4)));hold on;
                end
            end
            axes('Position',[0.55 0.4 0.4 0.55]);
            max_longi=(squeeze(max(img3,[],1)))';
            max_longi=[max_longi(:,end-hpsiz+1:end),max_longi,max_longi(:,1:hpsiz)];
            imagesc(max_longi,[-0.05,0.2]);colormap jet;hold on;
            if ~isempty(tracks)
                for itj=1:tracks(end,4)
                    trtmp=tracks(tracks(:,4)==itj,:);
                    xtr=trtmp(:,1);%-hpsiz;
                    ytr=trtmp(:,2)-hpsiz;
                    ttr=trtmp(:,3);
                    plot(xtr,ttr,'.-k');hold on;%,'Color',GenColor(itj/tracks(end,4)));hold on;
                end
            end
            xlabel('angle')
            ylabel('frames')
            set(gca,'XTick',(1:16:64)+hpsiz);
            set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
            title('longitutde max fluctuation kemograph')
            colorbar
            axes('Position',[0.05 0.05 0.9 0.3]);
            imagesc(std_img,[0.02 0.15]);colormap jet;axis image;colorbar;
            axis off;title('standard deviation map');hold on;
            if ~isempty(tracks)
                for itj=1:tracks(end,4)
                    trtmp=tracks(tracks(:,4)==itj,:);
                    xtr=trtmp(:,1)-hpsiz;
                    ytr=trtmp(:,2)-hpsiz;
                    plot(xtr,ytr,'.-k');hold on;%,'Color',GenColor(itj/tracks(end,4)));hold on;
                end
            end
%             pause
            print figure
            print(f,[path,'\pt_filter\',name,'_',num2str(inuc)],'-dpng');
        end
    end
end