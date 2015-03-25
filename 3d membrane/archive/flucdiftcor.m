% cell fluc and drift
clear all
close all
%set header
path='C:\nuclei\post analysis result_0.1';
datapath=[path,'\data'];
savepicdir=[path,'\pic_0.1'];

load([path,'\pt_filter\orientations.mat']);
load([path,'\pt_filter\alltracks.mat']);
nameNMBC={'wild_type','heh1','heh2','ima1','heh1heh2','heh1ima1','heh2ima1','heh1heh2ima1'};
nameMBC={'sp10_MBC','heh1_MBC','heh2_MBC','ima1_MBC','heh1heh2_MBC','heh1ima1_MBC','heh2ima1_MBC','heh1heh2ima1_MBC'};
nameAll=[nameNMBC,nameMBC];

% save all drift data
files=dir(datapath);
movienames={files(3:end).name}';

if ~exist([path,'\pt_filter\centerdrift.mat'],'file')
    nucdriftname=[];
    nucdrift=[];
    saveid=1;
    for ifile=3:length(files)
        display(['processing ',files(ifile).name]);
        load(fullfile(datapath,files(ifile).name));
        for inuc=1:nm.num_nuc
            centerdrift=zeros(nm.endframe,3);
            for iframe=1:nm.endframe
                centerdrift(iframe,:)=nm.nuclei{iframe,inuc}.origin_new;
            end
            nucdriftname{saveid}=[files(ifile).name(1:end-4),'_',num2str(inuc)];
            nucdrift{saveid}=centerdrift;
            saveid=saveid+1;
        end
    end
    save([path,'\pt_filter\centerdrift.mat'],'nucdriftname','nucdrift');
else
    load([path,'\pt_filter\centerdrift.mat']);
end

% parse index
nuctype=cell(size(nucind));
nucmovie=cell(size(nucind));
nucnumind=zeros(size(nucind));
for inuci=1:length(nucind)
    namefull=nucind{inuci};
    strtmpind=regexp(namefull,'_');
    nuctype{inuci}=namefull(1:strtmpind(end-1)-1);
    nucmovie{inuci}=namefull(1:strtmpind(end)-1);
    nucnumind(inuci)=str2double(namefull(strtmpind(end)+1:end));
end

%%
% extract locations
allpsi=cell(size(nameAll));
alltheta=cell(size(nameAll));
allcorr=cell(size(nameAll));
avgdriftstd=cell(size(nameAll));
driftspeed=cell(size(nameAll));
for iname=1:length(nameAll)
    name1=nameAll{iname};
    typeind=cellfun(@(x)(strcmp(x,name1)),nuctype);
    typetrj=trj(typeind);
    typemovie=nucmovie(typeind);
    typenucnumind=nucnumind(typeind);
    uniquemovie=unique(typemovie);
    psitmp=[];
    thetatmp=[];
    corrtmp1=[];
    driftstdtmp=[];
    driftspeedtmp=[];
    for imovie=1:length(uniquemovie)
        display(['processing ',uniquemovie{imovie}] );
        movieind=cellfun(@(x)(strcmp(x,uniquemovie{imovie})),typemovie);
        movietrj=typetrj(movieind);
        movienucnumind=typenucnumind(movieind);
        orientmp=allorientations{cellfun(@(x)strcmp(x(1:end-4),uniquemovie{imovie}),movienames)};
        drifttmp=nucdrift(cellfun(@(x)~isempty(strfind(x,uniquemovie{imovie})),nucdriftname));
        orientations=orientmp(movienucnumind);
        drifttmp2=drifttmp(movienucnumind);
        
        for inuctrj=1:length(movietrj)
            nuctrj=movietrj{inuctrj};
            centerdrift=drifttmp2{inuctrj};
            ncd=centerdrift-ones(size(centerdrift,1),1)*mean(centerdrift,1);
            otheta=orientations(inuctrj)/180*pi;
            ovec=[cos(otheta),sin(otheta)];
            ncdovec=(ncd(:,1)*ovec(1)+ncd(:,2)*ovec(2));
            driftstdtmp=[driftstdtmp;std(ncdovec)];
            dncdovec=[0;diff(ncdovec)];
            corrtmp2=[];
            for itrj=1:nuctrj(end,end)
                subtrj=nuctrj(nuctrj(:,end)==itrj,:);
                if ~isempty(subtrj)
                    if max(subtrj(:,4))>0.1
                        subtrjnonext=subtrj(subtrj(:,5)~=0,:);
                        %                         subtrjnonext=subtrj;
                        psis=(subtrjnonext(:,1)-8)/64*2*pi;
                        thetas=(subtrjnonext(:,2)-18)/64*2*pi;
                        heights=subtrjnonext(:,4);
                        ts=subtrjnonext(:,end-1);
                        driveovec=(heights.*cos(thetas).*cos(psis-otheta));
                        %                         crosscorr=@(x,y) mean(x.*y)/sqrt(mean(x.^2))/sqrt(mean(y.^2));
                        %                         cormatrix=crosscorr(dncdovec(ts),driveovec);
                        driftspeedtmp=[driftspeedtmp;abs(mean(dncdovec(ts)))];
                        cormatrix=corr(dncdovec(ts),driveovec);
                        corrtmp2=[corrtmp2;cormatrix(1)];
                        %                         plot(ts,driveovec*8*0.16,'r');hold on;
                        meanpsi=mean(psis);
                        meantheta=mean(thetas);
                        psi1= acos((cos(meanpsi-otheta))*cos(meantheta))/pi*180;
                        psitmp=[psitmp,psi1];
                        thetatmp=[thetatmp,meantheta];
                    end
                end
            end
            corrtmp1=[corrtmp1;corrtmp2];
            if strcmp(uniquemovie{imovie},'wild_type_04') &&  movienucnumind(inuctrj)==3;
                f1=figure;
                a1=axes();
                line((2:50)*2.5,dncdovec(2:50)*0.16*2.5,'color','b');
                ylabel('center velocity of nuclei (\mum/seconds)');
                set(a1,'XAxisLocation','bottom','YAxisLocation','left','YColor','b')
                set(a1,'Xtick',0:40:140);
                
                a2=axes('Position',get(a1,'Position'));
                line(ts*2.5,driveovec*8*0.16,'color','r');
                ylabel('fluctuation height (\mum)');
                set(a2,'XAxisLocation','bottom','YAxisLocation','right','YColor','r')
                set(a2,'Xtick',0:40:140);
                xlabel('time(seconds)');
                title('Nucleus motion driven by microtubule induced fluctuations');
                FigureFormat(f1);
            end
        end
    end
    %     psitmp(psitmp>90)=180-psitmp(psitmp>90);
    allpsi{iname}=psitmp;
    alltheta{iname}=thetatmp;
    allcorr{iname}=corrtmp1;
    avgdriftstd{iname}=driftstdtmp;
    driftspeed{iname}=driftspeedtmp;
end
% save([path,'\pt_filter\loc_orientation.mat'],'allpsi','alltheta');
%% correlation
xval=-0.95:0.1:0.95;
countscorr=zeros(length(xval),length(allcorr));
meancorr=zeros(size(allcorr));
secorr=zeros(size(allcorr));
meandriftstd=zeros(size(avgdriftstd));
sedriftstd=zeros(size(avgdriftstd));
meandriftspeed=zeros(size(driftspeed));
sedriftspeed=zeros(size(driftspeed));
for i=1:length(allcorr)
    corri=allcorr{i};
    ns=hist(corri,xval);
    countscorr(:,i)=ns;
    meancorr(i)=mean(allcorr{i});
    secorr(i)=std(allcorr{i})/sqrt(length(allcorr{i}));
    meandriftspeed(i)=mean((driftspeed{i}));
    sedriftspeed(i)=std((driftspeed{i}))/sqrt(length(driftspeed{i}));
end
% drift velocity bar
% barh(meandriftspeed(1:8));hold on;
% herrorbar(meandriftspeed(1:8),1:8,sedriftspeed(1:8),'k');
% set(gca,'YTick',1:8,'YTicklabel',nameAll(1:8));
% xlabel('correlation value');
% title('correlation between nuclei center motion and fluctuation motion');

% meancorrelation bar
% barh(meancorr(1:8));hold on;
% herrorbar(meancorr(1:8),1:8,secorr(1:8),'k');
% % axis([-1 1 .5 8.5])
% set(gca,'YTick',1:8,'YTicklabel',nameAll(1:8));
% xlabel('correlation value');
% title('correlation between nuclei center motion and fluctuation motion');

%scatter
% hsvcolor=hsv;
% for i=1:8
% scatter(meandriftspeed(i),meancorr(i),'fill','markerfacecolor',hsvcolor(i*8,:))
% hold on;
% end
% legend(strrep(nameAll(1:8),'_',' '));
% xlabel('centroid velocity (\mum/sec)');
% ylabel('correlation with microtubule induced fluctuations');
% FigureFormat(gcf);
%% plotting correlation of all
xval=-0.95:0.2:0.95;
countscorr=zeros(length(xval),length(allcorr));
meancorr=zeros(size(allcorr));
secorr=zeros(size(allcorr));
for i=1:length(allcorr)
    corri=allcorr{i};
    ns=hist(corri,xval);
    countscorr(:,i)=ns;
    meancorr(i)=mean(allcorr{i});
    secorr(i)=std(allcorr{i})/sqrt(length(allcorr{i}));
end

% wt
otherid=4;
bar(xval,[countscorr(:,1)/sum(countscorr(:,1)),...
    ...sum(countscorr(:,2:4),2)/sum(sum(countscorr(:,2:4))),...
    sum(countscorr(:,2:8),2)/sum(sum(countscorr(:,2:8))),]*100,'rb');
xlabel('correlation value');
ylabel('percentage (%)');
title('distribution of fluctuations');
box off
axis([-1 1 0 40])
legend('Wild Type','single knockouts')
FigureFormat(gcf)

    
    % all

for i=1:16
    clf
    bar(xval,countscorr(:,i)/sum(countscorr(:,i))*100);
    xlabel('correlation value');
    ylabel('percentage (%)');
    title('distribution of fluctuations');
    box off
    axis([-1 1 0 60])
    legend(nameAll{i})
    FigureFormat(gcf)
    pause
end
%% plotting cell orientation
psiallNMBC=[];
psiallMBC=[];
thetaall=[];
for iname=1:length(nameAll)/2
    psiallNMBC=[psiallNMBC,allpsi{iname}];
end
for iname=9:length(nameAll)
    psiallMBC=[psiallMBC,allpsi{iname}];
end

for iname=1:length(nameAll)
    thetaall=[thetaall,alltheta{iname}];
end

figure
xval=5:10:175;
dpNMBC=hist(psiallNMBC,xval);
dpMBC=hist(psiallMBC,xval);
%     bar(xval',[dpNMBC',dpMBC'],'Stacked');
bar(xval',dpNMBC,'b');hold on;
bar(xval',dpMBC,'r');
box off
set(gca,'Ytick',0:20:100);
set(gca,'Xtick',0:45:180);
set(gca,'XaxisLocation','bottom','YaxisLocation','left');
legend('non MBC','MBC');
xlabel('angles between cell orientation and fluctuation location (degree)');
ylabel('number of fluctuations detected')
title('locational distribution of fluctuations' )
FigureFormat(gcf);
print(gcf,[savepicdir,'\locall'],'-dpng');

%% figure, individual angle distribution
allpsi1=allpsi;
for i=1:16
    allpsi1{i}(allpsi1{i}>90)=180-allpsi1{i}(allpsi1{i}>90);
end
for i=1:8
    xval=5:10:85;
    clf
    dpNMBC=hist(allpsi1{i},xval);
    dpMBC=hist(allpsi1{i+8},xval);
    %     bar(xval',[dpNMBC',dpMBC'],'Stacked');
    bar(xval',dpNMBC,'b');hold on;
    bar(xval',dpMBC,'r');
    axis([0 90 0 20]);
    box off
    set(gca,'Ytick',0:5:20);
    set(gca,'Xtick',0:45:180);
    set(gca,'XaxisLocation','bottom','YaxisLocation','left');
    legend(nameAll{i},nameAll{i+8});
    xlabel('angles between cell orientation and fluctuation location (degree)');
    ylabel('number of fluctuations detected')
    title('locational distribution of fluctuations' )
    FigureFormat(gcf);
    pause
end

%% average nuclei drift along cell orientation