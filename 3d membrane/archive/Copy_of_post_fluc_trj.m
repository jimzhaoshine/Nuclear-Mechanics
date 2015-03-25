% make all the kemographs for max longtitude and equator line, save to
% fig, save all the trajectories

%use 0.05 drift corrected data
% path='C:\nuclei\post analysis result_0.05';
path='C:\nuclei\post analysis result_0.1';
% important parameters
% peak finding threshold
peak_th=0.05;
% tracking parameter
tparam.dim=2;
tparam.maxdisp=5;
tparam.mem=2;
tparam.good=10;
% track extension parameter
trjextwz=30;
if ~exist([path,'\pt_filter'],'dir')
    mkdir([path,'\pt_filter']);
end

files=dir([path,'\data']);
f=figure('Position',[100 100 1200 900]);
nucind=cell(10000,1);
trj=cell(10000,1);
savingindex=1;
for ifile=3:length(files);
    filename=files(ifile).name;
    [~,name,ext]=fileparts(filename);
    display(['processing ',name]);
    if strcmp(ext,'.mat') %&& ~isempty(strfind(name,'heh1heh2'))
        load(fullfile(path,'\data',filename));
        if ifile==3
            [ points2,faces2,psi,theta] = setup_img( );
        end
        %% go through all the nuclei
        for inuc=1:nm.num_nuc
            %image extension for gaussian filter
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
            meannuc=centralband(meannuc,nm,points2,size(psi));
            mean_img=meannuc.img;
            mean_img_ext=[mean_img(:,end-hpsiz+1:end),mean_img,mean_img(:,1:hpsiz)];
            mean_img_ext=[zeros(hpsiz,size(mean_img_ext,2));...
                mean_img_ext;zeros(hpsiz,size(mean_img_ext,2))];

            % construct 3d image series and filtered image;
            time_img=zeros(size(psi,1)+2*hpsiz,size(psi,2)+2*extsiz+2*hpsiz,nm.endframe);
            for iframe=1:nm.endframe
                nuc=nm.nuclei{iframe,inuc};
                img=nuc.img;
                dimg=(img-mean_img);%./mean_img;
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

            % find peaks and params for each image
            pos=[];
            for iframe=1:nm.endframe
                %filtered img
                ext_img_filtered=squeeze(time_img_filtered(:,:,iframe));
                ext_img=squeeze(time_img(:,:,iframe));
                %origin img
                % tracks must have peak value higher than 0.05
                pks=pkfnd(ext_img_filtered/mean_img_ext,peak_th,hpsiz);
                
                %filter out near cap peaks
                if ~isempty(pks)
%                     pks=pks(pks(:,2)>=hpsiz+3 & pks(:,2)<=hpsiz+14,:);
                    pks=pks(pks(:,2)>=hpsiz+5 & pks(:,2)<=hpsiz+15,:);
                end
                % getting specs for the peaks
                cnt=[];
                if ~isempty(pks)
                    for ipeak=1:size(pks,1)
                        wimg_filtered=ext_img_filtered(pks(ipeak,2)+(-hpsiz:hpsiz),pks(ipeak,1)+(-hpsiz:hpsiz));
                        wimg=ext_img(pks(ipeak,2)+(-hpsiz:hpsiz),pks(ipeak,1)+(-hpsiz:hpsiz));
                        bw=wimg_filtered>0.5*wimg_filtered(hpsiz+1,hpsiz+1);
                        wimg_weight=zeros(size(wimg_filtered));
                        wimg_weight(bw)=wimg_filtered(bw);
                        [imy,imx]=meshgrid(pks(ipeak,2)+(-hpsiz:hpsiz),pks(ipeak,1)+(-hpsiz:hpsiz));
                        my=mean(mean(wimg_weight.*imy))/mean(wimg_weight(:));
                        mx=mean(mean(wimg_weight.*imx))/mean(wimg_weight(:));
                        cnt=[cnt;mx,my,wimg(hpsiz+1,hpsiz+1),wimg_weight(hpsiz+1,hpsiz+1),sqrt(sum(bw(:)))];
                    end
                    pos=[pos;cnt,iframe+zeros(size(cnt,1),1)];
                end
            end
            % connect the trajectories
            tracks=[];
            tracksfull=[];
            if ~isempty(pos)
                ytracks=yaotrack(pos,tparam);
                if ~isempty(ytracks)
                    tracks=ytracks(:,[1 2 end-1 end]);
                    tracksfull=ytracks;
                end
            end
            
            % extend all the trajectories
            tracks_ext=[];
            endframe=nm.endframe;
            if ~isempty(tracksfull)
                for itrj=1:tracksfull(end,end)
                    trjtmp=tracksfull(tracksfull(:,end)==itrj,:);
                    if ~isempty(trjtmp)
                    t_start=trjtmp(1,end-1);
                    t_end=trjtmp(end,end-1);
                    mx_start=trjtmp(1,1);
                    mx_end=trjtmp(end,1);
                    my_start=trjtmp(1,2);
                    my_end=trjtmp(end,2);
                    pre_ts=(max(1,t_start-trjextwz):t_start-1)';
                    post_ts=(t_end+1:min(endframe,t_end+trjextwz))';
                    pre_x=zeros(size(pre_ts))+round(mx_start);
                    pre_y=zeros(size(pre_ts))+round(my_start);
                    post_x=zeros(size(post_ts))+round(mx_end);
                    post_y=zeros(size(post_ts))+round(my_end);
                    pre_h=squeeze(time_img_filtered(round(my_start),round(mx_start),pre_ts));
                    post_h=squeeze(time_img_filtered(round(my_end),round(mx_end),post_ts));
                    pre_h0=squeeze(time_img(round(my_start),round(mx_start),pre_ts));
                    post_h0=squeeze(time_img(round(my_end),round(mx_end),post_ts));
                    pre_w=zeros(size(pre_ts));
                    post_w=zeros(size(post_ts));
                    pre_ind=zeros(size(pre_ts))+itrj;
                    post_ind=zeros(size(post_ts))+itrj;
                    tracks_ext=[tracks_ext;pre_x,pre_y,pre_h0,pre_h,pre_w,pre_ts,pre_ind;...
                        trjtmp;post_x,post_y,post_h0,post_h,post_w,post_ts,post_ind];
                    end
                end
            end
            
            %save all the data to trjactories
            nucind{savingindex}=[name,'_',num2str(inuc)];
            trj{savingindex}=tracks_ext; 
            savingindex=savingindex+1;

            % plot
            clf
            % centerline kemograph
            axes('Position',[0.05 0.4 0.4 0.55]);
            img_centerline=(squeeze(time_img_filtered(10,:,:)))';
%             img_centerline=[img_centerline(:,end-hpsiz+1:end),img_centerline,img_centerline(:,1:hpsiz)];
            imagesc(img_centerline,[-0.05,0.2]);colormap jet;
            xlabel('angle')
            ylabel('frames')
            set(gca,'XTick',(1:16:64)+hpsiz);
            set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
            title('longitutde centerline fluctuation kemograph')
            colorbar;hold on;
            if ~isempty(tracks_ext)
                for itj=1:tracks_ext(end,end)
                    trtmp=tracks_ext(tracks_ext(:,end)==itj,:);
                    if ~isempty(trtmp)
                        xtr=trtmp(:,1);%-hpsiz;
                        ytr=trtmp(:,2);%-hpsiz;
                        ttr=trtmp(:,end-1);
                        plot(xtr,ttr,'.-k');hold on;%,'Color',GenColor(itj/tracks(end,4)));hold on;
                        
                    end
                end
            end
            % max kemograph
            axes('Position',[0.55 0.4 0.4 0.55]);
            max_longi=(squeeze(max(time_img_filtered,[],1)))';
%             max_longi=[max_longi(:,end-hpsiz+1:end),max_longi,max_longi(:,1:hpsiz)];
            imagesc(max_longi,[-0.05,0.2]);colormap jet;hold on;
            xlabel('angle')
            ylabel('frames')
            set(gca,'XTick',(1:16:64)+hpsiz);
            set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
            title('longitutde max fluctuation kemograph')
            colorbar            
            if ~isempty(tracks)
                for itj=1:tracks_ext(end,end)
                    trtmp=tracks_ext(tracks_ext(:,end)==itj,:);
                    if ~isempty(trtmp)
                    xtr=trtmp(:,1);%-hpsiz;
                    ytr=trtmp(:,2);%-hpsiz;
                    ttr=trtmp(:,end-1);
                    plot(xtr,ttr,'.-k');hold on;%,'Color',GenColor(itj/tracks(end,4)));hold on;
                    end
                end
            end

            % std map
            axes('Position',[0.05 0.05 0.6 0.3]);
            imagesc(std_img,[0.02 0.15]);colormap jet;axis image;colorbar;
            axis off;title('standard deviation map');hold on;
            if ~isempty(tracks)
                 for itj=1:tracks_ext(end,end)
                    trtmp=tracks_ext(tracks_ext(:,end)==itj,:);
                    if ~isempty(trtmp)
                    xtr=trtmp(:,1);%-hpsiz;
                    ytr=trtmp(:,2);%-hpsiz;
                    plot(xtr,ytr,'.-k');hold on;%,'Color',GenColor(itj/tracks(end,4)));hold on;
                    end
                 end
            end
            % plot all the tracks
            axes('Position',[0.7 0.05 0.25 0.3]);
            
            if ~isempty(tracks_ext)
                for itrj=1:tracks_ext(end,end)
                    trjtmp=tracks_ext(tracks_ext(:,end)==itrj,:);
                    trjh=trjtmp(:,4);
                    trjh0=trjtmp(:,3);%unfiltered center is very noisy, filtered result is good for time analysis
                    plot(trjh);hold on;
                end
            end
            %             pause
            print figure
            print(f,[path,'\pt_filter\',name,'_',num2str(inuc)],'-dpng');
        end
     end
end
nucind=nucind(cellfun(@(x)~isempty(x),trj));
trj=trj(cellfun(@(x)~isempty(x),trj));

save([path,'\pt_filter\','alltracks.mat'],'trj','nucind')