%% setup params
clear all;clc;close all;
run setup_header5.m;
rootpath='C:\nuclei\post analysis result_0.2';

strain_struct=struct('name',[],...
    'nuclei',[],... saved nuclei
    'numnuc',[],... number of nuclei
    'numfluc',[],... number of fluctuations
    'numMTfluc',[],...
    'thermalfluc',[],...
    'flucs',[]... fluctuations
    );

nuclei_struct=struct('name',[],...
    'movieid',[],... nuclei type
    'nucid',[],... nuclei id in movie
    'centerpos',[],... center position (um)
    'centervelo',[],... center velocity (um/s)
    'orientation',[],... orientation of cell
    'celllength',[],... cell length (um)
    'cellwidth',[],... cell width (um)
    'cellvolume',[],... cell volume (um^3)
    'nucleusvolume',[],... nuclear volume (im^3)
    'flucs',[],... fluctuations
    'size',[],... radius of the nuclei (um)
    'meanimg',[],... mean img (um)
    'stdimg',[],... standard deviation of movie (um)
    ...'thermalbw',[],... binary img identify nonfluc locations
    'meanr',[],...mean 3d contour (um)
    'good',[]...
    );

fluc_struct=struct(...
    'time',[],... time (s)
    'longitude',[],... longitude (rad)
    'latitude',[],... latitude (rad)
    'height',[],... height (um)
    'height_relative',[],...
    'width',[],... width (rad)
    'confidentind',[],... confident tracking (frames)
    'maxheightind',[],... maxheight peak index (frames)
    'maxheight',[],... maximum height from baseline (um)
    'maxheighttime',[],...maximum height timepoint (s)
    'maxheightlong',[],... maximum height longitude (rad)
    'maxheightlat',[],... maximum height latitude (rad)
    'avglong',[],... average longitude (rad)
    'avglat',[],... average latitude (rad)
    'meanwidth',[],... mean width (rad)
    'risetime',[],... risetime (s)
    'falltime',[],... falltime (s)
    'baseline',[],... %fluctuation baseline (um);
    'noisydata',[],... %how noisy is rise/fall time fitting
    'size',[],... %nuclei size
    'good',[]...
    );

% load good data
load([rootpath,'\goodnuclei.mat']);
% movienames1={files(3:end).name}';
% nuctypes=cell(size(gnmovienames));
% movieids=zeros(size(gnmovienames));
% for i=1:length(gnmovienames);
%     nuctypes{i}=gnmovienames{i}(1:end-7);
%     movieids(i)=str2double(gnmovienames{i}(end-5:end-4));
% end
%
% grouping basics
points=TriSphere(3);
hpsiz=16;
for itype=1:length(nameAll)
    name=nameAll{itype};
    typeind=find(strcmp(gnmovienames(:,1),name));
    nucsave=1;
    clear nucleitmp;
    %     nucleitmp=[];
    for imovie=typeind'
%         movieid=movieids(imovie);
        movieid=gnmovienames{imovie,3};
        load(fullfile(rootpath,'data',gnmovienames{imovie,1},gnmovienames{imovie,3}));
        numnuc=nm.num_nuc;
        zxr=nm.vox/nm.pix*nm.aberation;
        display(['processing ',gnmovienames{imovie,3}]);
        for inuc=1:numnuc
            nuctmp=nuclei_struct;
            %moviename movieid
            nuctmp.name=name;
            nuctmp.movieid=movieid;
            nuctmp.nucid=inuc;
            %orientation
            theta=nm.orientation(inuc)/180*pi;
            nuctmp.celllength=nm.celllength(inuc)*p2um;
            nuctmp.cellwidth=nm.cellwidth(inuc)*p2um;
            nuctmp.orientation=theta;
            nuctmp.cellvolume=(nuctmp.celllength-nuctmp.cellwidth)*pi/4*nuctmp.cellwidth^2+...
                +pi/6*(nuctmp.cellwidth).^3;
            nuctmp.good=goodnuclei{imovie}(inuc);
            %           right rotation
            rot_mat=[cos(theta), sin(theta),  0;...
                -sin(theta), cos(theta),  0;...
                0,           0,       1];
            %mean img, mean r, centerdift,
            sum_r=0;
%             centerdrift=zeros(nm.endframe,3);
            centerdriftraw=zeros(nm.endframe,3);
            nucsize=zeros(nm.endframe,1);
            [s1,s2]=size(nm.nuclei{1}.img);
            timeimg=zeros(s1,s2,nm.endframe);
            radii=zeros(nm.endframe,size(points,1));
            for iframe=1:nm.endframe
                nuc=nm.nuclei{iframe,inuc};
                %mean contour
                sum_r=sum_r+nuc.r_new;
                % time img
                timeimg(:,:,iframe)=nuc.img;
                %center drift
                cnt=nuc.origin_new;
%                 centerdrift(iframe,:)=(rot_mat*cnt')';
                centerdriftraw(iframe,:)=[nuc.origin_new(1:2),nuc.origin_new(3)*zxr];
                %nuclei radius
                radii(iframe,:)=nuc.r_new';
                % nuclei size
                if cnt(3)>=10
                    nucsize(iframe)=nuc.contour(10).area;
                elseif cnt(3)<=1
                    nucsize(iframe)=nuc.contour(1).area;
                else
                    lowerstack=floor(cnt(3));
                    upperstack=floor(cnt(3))+1;
                    upperweight=cnt(3)-lowerstack;
                    lowerweight=1-upperweight;
                    lowerarea=nm.nuclei{iframe,inuc}.contour(lowerstack).area;
                    upperarea=nm.nuclei{iframe,inuc}.contour(upperstack).area;
                    nucsize(iframe)=sqrt((lowerarea*lowerweight+upperarea*upperweight)/pi);
                end
            end
            mean_r=sum_r/nm.endframe;
            meannuc=nm.nuclei{1,inuc};
            meannuc.r_new=mean_r;
            meannuc=centralband(meannuc,nm,points2,size(psi),1);
            stdimg=std(timeimg,1,3);
            
            nuctmp.meanimg=meannuc.img*p2um;
            nuctmp.meanr=mean_r*p2um;
%             nuctmp.centerpos=centerdrift*p2um;
            nuctmp.centerposraw=centerdriftraw*p2um;
%             nuctmp.centervelo=[0 0 0;diff(centerdrift,1,1)]*p2um/f2s;
            nuctmp.size=mean(nucsize)*p2um;
            nuctmp.nucleusvolume=4*pi/3*nuctmp.size.^3;
            nuctmp.stdimg=stdimg*p2um;
            nuctmp.rmsf=std(radii)*p2um;
            
            % flucutations
            nucfluc=nm.trj{inuc};
            if ~isempty(nucfluc)
                numtrj=nucfluc(end,end);
                clear fluc;
                ifluc=1;
                for itrj=1:nucfluc(end,end)
                    subtrj=nucfluc(nucfluc(:,end)==itrj,1:end);
                    if ~isempty(subtrj)
                        subfluc=fluc_struct;
                        subfluc.time=subtrj(:,end-1)*f2s;
%                         subfluc.longitude=(subtrj(:,1)-hpsiz)/64*2*pi-theta;
%                         subfluc.latitude=(subtrj(:,2)-hpsiz-10)/64*2*pi;
                        subfluc.longitude=(subtrj(:,1)-hpsiz)/64*2*pi;
                        subfluc.latitude=(subtrj(:,2)-hpsiz-10)/64*2*pi;
                        subfluc.height=subtrj(:,4)*p2um;
                        subfluc.width=subtrj(:,5)/64*2*pi;
                        subfluc.meanwidth=mean(subfluc.width);
                        subfluc.height_relative=subfluc.height/nuctmp.size;
                        subfluc.size=nuctmp.size;
                        subfluc.confidentind=find(subtrj(:,5)~=0);
                        subfluc.avglong=mod(mean(subfluc.longitude(subfluc.confidentind)),2*pi);
                        subfluc.avglat=mean(subfluc.latitude(subfluc.confidentind));
                        subfluc.meanheight=mean(subfluc.height(subfluc.confidentind));
                        subfluc.good=nuctmp.good;
                        subfluc.inplane=abs(subfluc.avglat)<30/180*pi;
                        
                        influencezone=zeros(size(points(:,1)));
                        for ipts=1:length(subfluc.longitude)
                            psitmp=subfluc.longitude(ipts)+nuctmp.orientation;
                            thetatmp=subfluc.latitude(ipts);
                            cosdist=points(:,1)*cos(psitmp)*cos(thetatmp)...
                                +points(:,2)*sin(psitmp)*cos(thetatmp)...
                                +points(:,3)*sin(thetatmp);
                            influencezone(cosdist>cos(subfluc.width(ipts))*1.5)=1;
                        end
                        subfluc.influencezone=influencezone;
                        
                        fluc(ifluc)=subfluc;
                        ifluc=ifluc+1;
                    end
                end
                nuctmp.flucs=fluc;
            else
                nuctmp.flucs=[];
            end
            
            %get non fluc zone for this nuclei
            fluczone=zeros(size(points(:,1)));
            for ifluc=1:length(nuctmp.flucs)
                fluczone=fluczone+nuctmp.flucs(ifluc).influencezone;
            end
            nuctmp.nonfluczone= fluczone==0;
            %             if max(centerdrift(:,3))<=8 && min(centerdrift(:,3))>=3
            nucleitmp(nucsave)=nuctmp;
            nucsave=nucsave+1;
            %             end
        end
    end 
    %         rmsfall=[];
    %         zrange=find(abs(points(:,3))<0.5);
    %         for jj=1:length(nucleitmp)
    %             rmsfall=[rmsfall,nucleitmp(jj).rmsf(zrange)];
    %         end
    %            bins=0:0.0001:0.15;
    %         [counts]=hist(rmsfall,bins);
    %         cumcounts=cumsum(counts)./sum(counts);
    %         plot(bins,cumcounts,'linewidth',2);hold on;
    %         xlabel('rmsf(\mum)');ylabel('cummulative probability');axis([0 0.15 0 1])
    
        %save to strains
    strainall(itype)=strain_struct;
    strainall(itype).name=name;
    if exist('nucleitmp')
        if ~isempty(nucleitmp)
            strainall(itype).nuclei=nucleitmp;
            strainall(itype).numnuc=length(nucleitmp);
        else
            warning(['no nuceli available for ',name]);
        end
    end
end

% time analysis
for itype=1:length(nameAll)
    display(['processing ',nameAll{itype}]);
    name=strainall(itype).name;
    nuclei=strainall(itype).nuclei;
    for inuc=1:length(nuclei)
        flucs=nuclei(inuc).flucs;
        for ifluc=1:length(flucs)
            subfluc=flucs(ifluc);
            ts=subfluc.time;
            hs=subfluc.height;
            [mh,mhind]=max(subfluc.height);
            fitfun=@(P)risefall(P,ts)-hs;
            p0=[ts(mhind),50,50,0.3,-0.1];
            lb=[ts(1),10,10,0.05,-0.15];
            ub=[ts(end),min(ts(end)-ts(1),150),...
                min(ts(end)-ts(1),150),0.5,0.05];
            options=optimset('Display','off');
            p=lsqnonlin(fitfun,p0,lb,ub,options);
            
            if p(1)-ts(1)<=15
                noisydata='rise';
            elseif ts(end)-p(1)<=15
                noisydata='fall';
            else
                noisydata='none';
            end
            subfluc.noisydata=noisydata;
            subfluc.maxheighttime=p(1);
            subfluc.risetime=p(2);
            subfluc.falltime=p(3);
            subfluc.maxheight=p(4);
            subfluc.baseline=p(5);
            [~,subfluc.maxheightind]=min(abs(p(1)-ts));
            subfluc.maxheightlong=mod((subfluc.longitude(subfluc.maxheightind)),2*pi);
            subfluc.maxheightlat=mod((subfluc.latitude(subfluc.maxheightind)),2*pi);
            strainall(itype).nuclei(inuc).flucs(ifluc)=subfluc;
            %             clf
            %             plot(ts,hs);hold on;
            %             plot(ts,risefall(p,ts),'r');hold on;
            %             title(subfluc.noisydata);
            %             axis([0 250 -0.15 0.3])
            %             text(200,0.2,{['risetime ',num2str(subfluc.risetime)],...
            %                 ['falltime ',num2str(subfluc.falltime)]})
            %             FigureFormat(gcf);
            %             pause
        end
    end
end
% center motion and centroid correlation
% for itype=1
%     for inuc=1:length(strainall(itype).nuclei)
%         nuc=strainall(itype).nuclei(inuc);
%         for ifluc=1:length(nuc.flucs)
%             fluc=nuc.flucs(ifluc);
%             cnt=nuc.centerpos-ones(101,1)*mean(nuc.centerpos,1);
%             hx=fluc.height*cos
%             clf
%             plot((1:101)*2.5,cnt);legend('x','y','z');hold on;
%             plot(fluc.time,'b');hold on;
%             pause
%         end
%     end
% end

% resave flucs
for itype=1:length(nameAll)
    nuclei=strainall(itype).nuclei;
    flucs=[];
    for inuc=1:length(nuclei)
        flucs=[flucs,nuclei(inuc).flucs];
    end
    strainall(itype).flucs=flucs;
    strainall(itype).numfluc=length(flucs);
end
save([rootpath,'\strainall.mat'],'strainall');
load([rootpath,'\strainall.mat']);
