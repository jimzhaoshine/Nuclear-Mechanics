%% setup params
clear all;
clc;
close all;
rpath='C:\nuclei\post analysis result_0.2';
nameNMBC={'wild_type','heh1','heh2','ima1','heh1heh2','heh1ima1','heh2ima1','heh1heh2ima1'};
nameMBC={'sp10_MBC','heh1_MBC','heh2_MBC','ima1_MBC','heh1heh2_MBC','heh1ima1_MBC','heh2ima1_MBC','heh1heh2ima1_MBC'};
nameIntercleave=reshape([nameNMBC;nameMBC],1,16);
nameAll=[nameNMBC,nameMBC];
colorAll=[255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162;255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162]/255;

numM=length(nameNMBC);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
sd=[rpath,'\result'];
mkdir(sd)

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

% % load all the orientations
% files=dir([path1,'\data']);
% load('C:\nuclei\data\loc.mat');
% % movienames1={files(3:end).name}';
% nuctypes=cell(size(movienames1));
% movieids=zeros(size(movienames1));
% for i=1:length(movienames1);
%     nuctypes{i}=movienames1{i}(1:end-7);
%     movieids(i)=str2double(movienames1{i}(end-5:end-4));
% end

% load all trjectories
load([rpath,'\pt_filter\alltracks.mat']);
fluctypes=cell(size(nucind));
flucmovieids=zeros(size(nucind));
flucinucs=zeros(size(nucind));
for i=1:length(nucind)
    splitter=regexp(nucind{i},'_');
    fluctypes{i}=nucind{i}(1:splitter(end-1)-1);
    flucmovieids(i)=str2double(nucind{i}(splitter(end-1)+1:splitter(end)-1));
    flucinucs(i)=str2double(nucind{i}(splitter(end)+1:end));
end

% load good data
load([rpath,'\goodnuclei.mat']);
% movienames1={files(3:end).name}';
nuctypes=cell(size(gnmovienames));
movieids=zeros(size(gnmovienames));
for i=1:length(gnmovienames);
    nuctypes{i}=gnmovienames{i}(1:end-7);
    movieids(i)=str2double(gnmovienames{i}(end-5:end-4));
end

%% grouping basics
points=TriSphere(3);
hpsiz=16;
for itype=1:length(nameAll)
    name=nameAll{itype};
    typeind=find(strcmp(nuctypes,name));
    nucsave=1;
    clear nucleitmp;
    for imovie=typeind'
        movieid=movieids(imovie);
        load([rpath,'\data\',gnmovienames{imovie}]);
        numnuc=nm.num_nuc;
        display(['processing ',gnmovienames{imovie}]);
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
            centerdrift=zeros(nm.endframe,3);
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
                centerdrift(iframe,:)=(rot_mat*cnt')';
                centerdriftraw(iframe,:)=nuc.center;
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
            nuctmp.centerpos=centerdrift*p2um;
            nuctmp.centerposraw=centerdriftraw*p2um;
            nuctmp.centervelo=[0 0 0;diff(centerdrift,1,1)]*p2um/f2s;
            nuctmp.size=mean(nucsize)*p2um;
            nuctmp.nucleusvolume=4*pi/3*nuctmp.size.^3;
            nuctmp.stdimg=stdimg*p2um;
            nuctmp.rmsf=std(radii)*p2um;
            
            % flucutations
            findtrj=find(strcmp(fluctypes,name)&flucmovieids==movieid&flucinucs==inuc);
            nucfluc=trj{findtrj};
            if ~isempty(nucfluc)
                numtrj=nucfluc(end,end);
                clear fluc;
                ifluc=1;
                for itrj=1:nucfluc(end,end)
                    subtrj=nucfluc(nucfluc(:,end)==itrj,1:end);
                    if ~isempty(subtrj)
                        subfluc=fluc_struct;
                        subfluc.time=subtrj(:,end-1)*f2s;
                        subfluc.longitude=(subtrj(:,1)-hpsiz)/64*2*pi-theta;
                        subfluc.latitude=(subtrj(:,2)-hpsiz-10)/64*2*pi;
                        subfluc.height=subtrj(:,4)*p2um;
                        subfluc.width=subtrj(:,5)/64*2*pi;
                        subfluc.meanwidth=mean(subfluc.width);
                        subfluc.height_relative=subfluc.height/nuctmp.size;
                        subfluc.size=nuctmp.size;
                        subfluc.confidentind=find(subtrj(:,5)~=0);
                        subfluc.avglong=mod(mean(subfluc.longitude(subfluc.confidentind)),2*pi);
                        subfluc.avglat=mod(mean(subfluc.latitude(subfluc.confidentind)),2*pi);
                        subfluc.good=nuctmp.good;
                        fluc(ifluc)=subfluc;
                        ifluc=ifluc+1;
                    end
                end
                nuctmp.flucs=fluc;
            else
                nuctmp.flucs=[];
            end
%             if max(centerdrift(:,3))<=8 && min(centerdrift(:,3))>=3
                nucleitmp(nucsave)=nuctmp;
                nucsave=nucsave+1;
%             end
        end
    end
    
    %save to strains
    strainall(itype)=strain_struct;
    strainall(itype).name=name;
    if ~isempty(nucleitmp)
        strainall(itype).nuclei=nucleitmp;
        strainall(itype).numnuc=length(nucleitmp);
    else
        warning(['no nuceli available for ',name]);
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
            
%             
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

% centroid correlation
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
save([rpath,'\pt_filter\strainall.mat'],'strainall');
load([rpath,'\pt_filter\strainall.mat']);
%% rmsf plot
zrange=find(abs(points(:,3))<=0.5);
xval=0.01:0.005:0.1;
rmsfcounts=zeros(length(strainall),length(xval));
for itype=1:length(strainall)
    str=strainall(itype);
    rmsftmp=[];
    for inuc=1:length(str.nuclei)
        rmsftmp=[rmsftmp,str.nuclei(inuc).rmsf];
    end
    rmsfctmp=hist(rmsftmp,xval);
    rmsfcounts(itype,:)=rmsfctmp/sum(rmsfctmp);
end

rmsfchose=1:16;%[1 9];
plot(xval'*ones(size(rmsfchose)), rmsfcounts(rmsfchose,:)');
legend(nameAll(rmsfchose));

%% nuclei radius size
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp=[strainall(itype).nuclei.size];
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
title('Nuclear Size');
axis([1 1.4 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('radius(\mum)');
FigureFormat(gcf)
print(gcf,[sd,'\size'],'-dpng');
%% cell length width ratio
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp1=[strainall(itype).nuclei.celllength];
    datatmp2=[strainall(itype).nuclei.cellwidth];
    datatmp=datatmp1./datatmp2;
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
title('Cell length - cell width ratio');
axis([1 4 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('ratio');
FigureFormat(gcf)
print(gcf,[sd,'\lengthwidthratio'],'-dpng');
%% nuclei volume
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp=[strainall(itype).nuclei.nucleusvolume];
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
title('Nuclear Volume');
axis([1 15 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('volume (\mum^3)');
FigureFormat(gcf)
print(gcf,[sd,'\nucvol'],'-dpng');
%% nucleus cell ratio
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp1=[strainall(itype).nuclei.cellvolume];
    datatmp2=[strainall(itype).nuclei.nucleusvolume];
    datatmp=datatmp2./datatmp1;
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
title('Nucleus to Cell ratio');
axis([0 0.12 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('N/C ratio');
FigureFormat(gcf)
print(gcf,[sd,'\nucleuscellratio'],'-dpng');

%% cell volume
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp=[strainall(itype).nuclei.cellvolume];
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
title('Cell Volume');
% axis([0 7200 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('cell volume(\mum^3)');
FigureFormat(gcf)
print(gcf,[sd,'\cellvol'],'-dpng');

%% size, length correlation
colortmp=jet(4);
for itype=1:4
    datatmp1=[strainall(itype).nuclei.nucleusvolume];
    datatmp2=[strainall(itype).nuclei.cellvolume];
    plot(datatmp2,datatmp1,'o','Color',colortmp(itype,:));hold on;
end
ylabel('nuclear volumes (\mum^3)');
xlabel('cell volumes (\mum^3)');
% axis([ 2 2 17]);
legend(strrep(nameAll(1:4),'_',' '),'Location','Northwest');
FigureFormat(gcf)
print(gcf,[sd,'\nucvolvsvol'],'-dpng');
%% fit each nuclei
datalength=16;
fitval=zeros(datalength,3);
for itype=1:datalength
    datatmp1=[strainall(itype).nuclei.nucleusvolume];
    datatmp2=[strainall(itype).nuclei.cellvolume];
    P=lsqnonlin(@(P)P*datatmp2-datatmp1,0.1);
%     [p,S]=polyfit(datatmp2,datatmp1,1);
    p=[P,0];
    yresid=datatmp1-polyval(p,datatmp2);
    SSresid = sum(yresid.^2);
    SStotal = (length(datatmp1)-1) * var(datatmp1);
    rsq=1-SSresid/SStotal;
    rsq
    fitval(itype,:)=[p,sqrt(rsq)];
    %clf
    plot(datatmp2,datatmp1,'.','Color',colorAll(itype,:));hold on;
    plot(0:200,polyval(p,0:200),'-k')
    xlabel('cell volume (\mum^3)');
    ylabel('nuclei volume (\mum^3)');
    axis([ 0 200 0 20]);
    legend(strrep(nameAll(itype),'_',' '),'Location','Northwest');
    FigureFormat(gcf)
%     pause
%     barh(datalength+1-itype,p(1),'facecolor',colorAll(itype,:));hold on;
end
% title('Base cell size');
% axis([0 2 0 datalength+1]);
% set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
% xlabel('volume(\mum^3)');
% FigureFormat(gcf)

%% cell length?
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp=[strainall(itype).nuclei.celllength];
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
title('Cell Length');
axis([6 12 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('radius(\mum)');
FigureFormat(gcf)
print(gcf,[sd,'\celllength'],'-dpng');
%% cell width?
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1001)
for itype=1:datalength
    datatmp=[strainall(itype).nuclei.cellwidth];
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
title('Cell Width');
axis([1 5 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('radius(\mum)');
FigureFormat(gcf)
print(gcf,[sd,'\cellwidth'],'-dpng');
%% nuclei center drift
datalength=16;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
figure(1011)
for itype=1:datalength
    cnts={strainall(itype).nuclei.centerpos};
    datatmpx=cellfun(@(x)x(:,1),cnts,'UniformOutput',0);
    datatmp=cellfun(@(x)std(polyval(polyfit(1:101,x',1),1:101)-x'),datatmpx);
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
    set(get(B,'child'),'facea',.1)
end
herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;

for itype=1:datalength
    cnts={strainall(itype).nuclei.centerpos};
    datatmpy=cellfun(@(x)x(:,2),cnts,'UniformOutput',0);
    datatmp=cellfun(@(x)std(polyval(polyfit(1:101,x',1),1:101)-x'),datatmpy);
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
    set(get(B,'child'),'facea',.5)

end
herrorbar(meandata,datalength:-1:1,sedata,'w');hold on;

for itype=1:datalength
    cnts={strainall(itype).nuclei.centerpos};
    datatmpz=cellfun(@(x)x(:,3),cnts,'UniformOutput',0);
    datatmp=cellfun(@(x)std(polyval(polyfit(1:101,x',1),1:101)-x'),datatmpz);
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
    set(get(B,'child'),'facea',1)
end
herrorbar(meandata,datalength:-1:1,sedata,'w');hold on;

title('Nuclear Center Motion');
axis([0 .15 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('Standard Deviation of Position(\mum)');
FigureFormat(gcf)
print(gcf,[sd,'\centerdrift'],'-dpng');

%% fluctuations per nuc per min

meandata=zeros(16,1);
sedata=zeros(16,1);
stddata=zeros(16,1);
figure(1002)
fparam=[0 0];
for itype=1:length(nameAll)
    [ s_ind, ds_ind ] = SelectMtFluc( strainall(itype).flucs,fparam);
    datatmp=cellfun(@(x)sum(SelectMtFluc(x,fparam)),...
        {strainall(itype).nuclei.flucs})/0.5/250*60;
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    B=barh(17-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        set(get(B,'child'),'facea',.3)
end
herrorbar(meandata,16:-1:1,sedata,'k');hold on;
fparam=[.15 50];
for itype=1:length(nameAll)
    [ s_ind, ds_ind ] = SelectMtFluc( strainall(itype).flucs,fparam);
    datatmp=cellfun(@(x)sum(SelectMtFluc(x,fparam)),...
        {strainall(itype).nuclei.flucs})/0.5/250*60;
    meandata(itype)=mean(datatmp);
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    barh(17-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
herrorbar(meandata,16:-1:1,sedata,'w');hold on;
title(' Fluctuation Frequency');
axis([0 1.3 0 17]);
set(gca,'YTick',1:16,'YTicklabel',nameAll(16:-1:1));
xlabel('number per min per nuclei (min^{-1})');
FigureFormat(gcf)
print(gcf,[sd,'\flucpermin'],'-dpng');


%% risetime and falltime of fluctuations
datalength=8;
meandata=zeros(datalength,1);
sedata=zeros(datalength,1);
stddata=zeros(datalength,1);
f=figure(1003);
set(f,'Position',[50 50 800 600]);
fparam=[0.15 50];
for itype=1:datalength
    [ s_ind, ds_ind ] = SelectMtFluc( strainall(itype).flucs,fparam);
    datatmp=[strainall(itype).flucs.risetime];
    realrise=find(~strcmp({strainall(itype).flucs.noisydata},'rise') & s_ind);
    meandata(itype)=mean(datatmp(realrise));
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    B=barh(datalength+1-itype,meandata(itype),'facecolor',...
        colorAll(itype,:));hold on;
    set(get(B,'child'),'facea',.3)
end
h1=herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
for itype=1:datalength
    [ s_ind, ds_ind ] = SelectMtFluc( strainall(itype).flucs,fparam);
    datatmp=[strainall(itype).flucs.falltime];
    realfall=find(~strcmp({strainall(itype).flucs.noisydata},'fall')& s_ind);
    meandata(itype)=mean(datatmp(realfall));
    stddata(itype)=std(datatmp);
    sedata(itype)=std(datatmp)./sqrt(length(datatmp));
    B=barh(datalength+1-itype,meandata(itype),'facecolor',...
        colorAll(itype,:));hold on;
end
h2=herrorbar(meandata,(datalength:-1:1),sedata,'w');hold on;
lh=legend([h1,h2],'rise time','fall time');
title(' Rise time and fall time of Fluctuations');
axis([0 100 0 datalength+1]);
set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
xlabel('risetime (s)');
FigureFormat(gcf)
set(lh, 'Color', [.8 .8 .9])
pause(.2)
print(f,[sd,'\risefall'],'-dpng');

%% number of nuclei

meandata=zeros(16,1);
sedata=zeros(16,1);
stddata=zeros(16,1);
figure(1004)
for itype=1:length(nameAll)
    meandata(itype)=length(strainall(itype).nuclei);
    barh(17-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
end
title(' Nuclei Collected');
axis([0 120 0 17]);
set(gca,'YTick',1:16,'YTicklabel',nameAll(16:-1:1));
xlabel('number');
FigureFormat(gcf)
print(gcf,[sd,'\number'],'-dpng');


%% fluctuation angle distribution
% load([path,'\pt_filter\strainall.mat'],'strainall');
for itype=1:length(nameAll)
    f=figure(1005);
    set(f,'Position',[50 50 800 600]);
    clf
    name=strainall(itype).name;
    nuclei=strainall(itype).nuclei;
    flucs=strainall(itype).flucs;
    xval=5:10:175;
    thetas=abs((mod([flucs.avglong]+pi,2*pi)-pi))/pi*180;
    counts=hist(thetas,xval);
    counts=counts./sum(counts);
    a1=axes('Position',[0.2 0.2 .6 .6]);
    axis([0 180 0 0.3]);
    box off;axis off;
    bh=bar(xval, counts,'r');hold on;
    set(a1,'Xtick',0:45:180,'visible','off');
    set(get(bh,'children'),'facea',.2);
    
    a2=axes('Position',[0.2 0.2 .6 .6]);
    scatter(thetas,[flucs.maxheight],'linewidth',2);hold on;
    axis([0 180 0 0.4]);
    set(a2,'Xtick',0:45:180,'Ytick',0:0.1:0.4,...
        'YaxisLocation','left','YColor','b');box off;
    ylabel('fluctuation max height (um)');
    
    a3=axes('Position',[0.2 0.2 .6 .6]);
    confidentdata=find(strcmp({strainall(itype).flucs.noisydata},'none'));
    scatter(thetas(confidentdata),([flucs(confidentdata).risetime]...
        +[flucs(confidentdata).falltime]),'linewidth',2,...
        'marker','*','markeredgecolor','g');
    axis([0 180 0 150]);
    
    set(a3,'Xtick',0:45:180,'Ytick',0:50:100,...
        'YaxisLocation','right','YColor','green');
    ylabel('fluctuation duration (s)');
    
    xlabel('angle between flucatuation and cell orientation (degree)');
    title(strrep(name,'_',' '));
    FigureFormat(gcf)
    print(gcf,[sd,'\fluc',name],'-dpng');
end
% try to differentiate mt from non mt and thermal,by height, location, time



%% scatter of duration and height
for itype=1:length(nameAll)
    datatmp1=[strainall(itype).flucs.risetime]+[strainall(itype).flucs.risetime];
    datatmp2=[strainall(itype).flucs.maxheight];
%     confidentdata=find(strcmp({strainall(itype).flucs.noisydata},'none'));
    scatter(datatmp2,datatmp1);
    axis([0.05 0.4 0 300])
    title(strrep(nameAll{itype},'_',' '))
    FigureFormat(gcf)
    pause

end

