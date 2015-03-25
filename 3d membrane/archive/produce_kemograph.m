% make all the kemographs for max longtitude and equator line, save to
% fig
clear all;
close all;
path='C:\nuclei\post analysis result_0.1';
files=dir([path,'\data']);
if ~exist([path,'\kemo'],'dir')
    mkdir([path,'\kemo']);
end
f=figure('Position',[100 100 1600 800]);
for i=3:length(files);
    filename=files(i).name;
    [~,name,ext]=fileparts(filename);
    display(name);
    if strcmp(ext,'.mat')
        load(fullfile([path,'\data'],filename));
        if i==3
            [ points2,faces2,psi,theta] = setup_img( );
        end
        for inuc=1:nm.num_nuc
            %% calculate mean contour
            sum_r=0;
            for iframe=1:nm.endframe
                sum_r=sum_r+nm.nuclei{iframe,inuc}.r_new;
            end
            mean_r=sum_r/nm.endframe;
            mean_pts=[mean_r.*nm.points(:,1),mean_r.*nm.points(:,2),mean_r.*nm.points(:,3)];
            meannuc=nm.nuclei{1,inuc};
            meannuc.r_new=mean_r;
            meannuc=centralband(meannuc,nm,points2,size(psi));
            mean_img=meannuc.img;
            time_img=zeros(size(psi,1),size(psi,2),nm.endframe);
            abs_time_img=zeros(size(psi,1),size(psi,2),nm.endframe);
            for iframe=1:nm.endframe
                nuc=nm.nuclei{iframe,inuc};
                img=nuc.img;
                dimg=img-mean_img;
                time_img(:,:,iframe)=dimg./mean_img;
            end
            % show kemograph
            clf
            subplot(1,2,1)
            % max fluc at diff longitude at time series
            imagesc(squeeze(time_img(10,:,:))',[-0.05,0.2]);colormap jet;
            xlabel('angle')
            ylabel('frames')
            set(gca,'XTick',10:16:74);
            set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
            title('longitutde centerline fluctuation kemograph')
            colorbar
            subplot(1,2,2)
            max_longi=squeeze(max(time_img,[],1));
            imagesc(max_longi',[-0.05,0.2]);colormap jet;
            xlabel('angle')
            ylabel('frames')
            set(gca,'XTick',10:16:74);
            set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
            title('longitutde max fluctuation kemograph')
            colorbar
            print(f,[path,'\kemo\',name,'_',num2str(inuc)],'-dpng');
        end
    end
end