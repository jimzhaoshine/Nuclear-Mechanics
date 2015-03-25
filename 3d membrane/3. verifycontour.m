%% verify movie quality
clear all;clc;close all;
rootpath='C:\nuclei\post analysis result_0.2';
verifypath=[rootpath,'\verification'];
mkdir(verifypath);
run setup_header3.m;

[points,faces,edges,neighbors]=TriSphere(3);
neighbors(1:12,6)=(1:12)';
zrange=find(abs(points(:,3))<0.5);
gnmovienames={};
goodnuclei=[];

allfiles=dir(fullfile(rootpath,'data'));
allnames={allfiles(3:end).name};
allcolors=jet(length(allnames));
% files=dir([rootpath,'\data\*.mat']);
% names={files.name};
% nuctypes=cellfun(@(x)x(1:end-7),names,'UniformOutput',0);
% movieids=cellfun(@(x)str2double(x(end-5:end-4)),names);
f1=figure(3001);
f2=figure(3002);
f3=figure(3003);
f4=figure(3004);
f5=figure(3005);
f6=figure(3006);
f7=figure(3007);
f8=figure(3008);
for itype=1:length(allnames)
    clf(f1);
    clf(f2);
    clf(f3);
    clf(f4);
    clf(f5);
    clf(f7);
    clf(f8);
%     typeind=find(strcmp(names,nameAll{itype}));
%     nucsave=1;
%     nummovie=length(typeind);
%     moviecolor=hsv(nummovie);
    typermsf=[];
%     legendids=cellfun(@num2str,num2cell(movieids(typeind)),'UniformOutput',0);
    moviefiles=dir(fullfile(rootpath,'data',allnames{itype},'*.mat'));
    movienames={moviefiles.name};
    legendids=movienames;
    moviecolor=jet(length(movienames));
    for imovie=1:length(movienames)
%         movieid=movieids(imovie);
%         load([rootpath,'\data\',names{imovie}]);
        load(fullfile(rootpath,'data',allnames{itype},movienames{imovie}));
        display(['processing ',movienames{imovie}]);
        goodnucleitmp=zeros(1,nm.num_nuc);
        %%
        rmsf=zeros(length(zrange),nm.num_nuc);
        dr2s=zeros(length(zrange),nm.num_nuc);
        drs=zeros(nm.num_nuc,length(zrange),nm.endframe);
        existflags=zeros(nm.num_nuc,nm.endframe);
        xs=zeros(nm.num_nuc,nm.endframe);
        ys=zeros(nm.num_nuc,nm.endframe);
        zs=zeros(nm.num_nuc,nm.endframe);
        ozs=zeros(nm.num_nuc,nm.endframe);
        dcs=zeros(nm.num_nuc,nm.endframe);
        for inuc=1:nm.num_nuc
            r_s=zeros(length(zrange),nm.endframe);
            dr_s=zeros(length(zrange),nm.endframe);
            for iframe=1:nm.endframe
                nuc=nm.nuclei{iframe,inuc};
                allr=nuc.r_new;
                neighbor_r=allr(neighbors);
                dr2=sum((allr*ones(1,6)-neighbor_r).^2,2)/6;

                r=nuc.r_new(zrange);
                r_s(:,iframe)=r;
                dr_s(:,iframe)=dr2(zrange);
                existflags(inuc,iframe)=nuc.exitflag;
                xs(inuc,iframe)=nuc.origin_new(1);
                ys(inuc,iframe)=nuc.origin_new(2);
                zs(inuc,iframe)=nuc.origin_new(3);
            end
%             drs(inuc,:,:)=r_s-mean(r_s,2)*ones(1,size(r_s,2));
            xs(inuc,:)=xs(inuc,:)-mean(xs(inuc,:));
            ys(inuc,:)=ys(inuc,:)-mean(ys(inuc,:));
            ozs(inuc,:)=zs(inuc,:);
            zs(inuc,:)=zs(inuc,:)-mean(zs(inuc,:));
            dcs(inuc,:)=sqrt(xs(inuc,:).^2+ys(inuc,:).^2+zs(inuc,:).^2)*p2um;
            rmsf(:,inuc)=std(r_s,1,2)*p2um;
            dr2s(:,inuc)=max(dr_s,[],2)';
            
            if max(ozs(inuc,:))<=8 && min(ozs(inuc,:))>=3 ...
                    && max(rmsf(:,inuc))<0.3 && mean(rmsf(:,inuc))<0.1 ...
                     && max(dcs(inuc,:))<0.6 ...
                     && max(dr2s(:,inuc))<0.5
                goodnucleitmp(inuc)=1;
            else
                goodnucleitmp(inuc)=0;
            end
        end
        gnmovienames=[gnmovienames;{allnames{itype},' ',movienames{imovie}}];
        goodnuclei=[goodnuclei;{goodnucleitmp}];
        
        figure(f1);
        plot(existflags(find(goodnucleitmp),:)','color',moviecolor(imovie,:));hold on;
        plot(existflags(find(~goodnucleitmp),:)','color',moviecolor(imovie,:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('existflag');axis([0 101 -2 4])

        figure(f2);
        plot(ozs(find(goodnucleitmp),:)','color',moviecolor((imovie),:));hold on;
        plot(ozs(find(~goodnucleitmp),:)','color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('z position (slice)');axis([0 101 0 11])
        
        figure(f3);
        plot(dcs(find(goodnucleitmp),:)','color',moviecolor((imovie),:));hold on;
        plot(dcs(find(~goodnucleitmp),:)','color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('drift (\mum)');axis([0 101 0 2])
        
        figure(f4)
        plot(rmsf(:,find(goodnucleitmp)),'color',moviecolor((imovie),:));hold on;
        plot(rmsf(:,find(~goodnucleitmp)),'color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('angles');ylabel('rmsf (\mum)');axis([0 length(zrange) 0 0.4])
        
        figure(f5)
        bins=0:0.0025:0.15;
        rmsf1=rmsf(:,find(goodnucleitmp));
        [counts]=hist(rmsf1(:),bins);
        cumcounts=cumsum(counts)./sum(counts);
        plot(bins,cumcounts,'linewidth',2,'color',moviecolor((imovie),:));hold on;
        xlabel('rmsf(\mum)');ylabel('cummulative probability');axis([0 0.15 0 1])
        legend(legendids);
        
        figure(f7)
        plot(dr2s(:,find(goodnucleitmp)),'color',moviecolor((imovie),:));hold on;
        plot(dr2s(:,find(~goodnucleitmp)),'color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('sum dr square over 6 (pixel^2)');axis([0 length(zrange) 0 1])
        
        figure(f8);set(f8,'Position',[0 0 1500 1000])
        subplot(2,3,1)
        plot(existflags(find(goodnucleitmp),:)','color',moviecolor((imovie),:));hold on;
        plot(existflags(find(~goodnucleitmp),:)','color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('existflag');axis([0 101 -2 4]);title('exitflag');
        subplot(2,3,2)
        plot(ozs(find(goodnucleitmp),:)','color',moviecolor((imovie),:));hold on;
        plot(ozs(find(~goodnucleitmp),:)','color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('z position (slice)');axis([0 101 0 11]);title('zposition');
        subplot(2,3,3)
        plot(dcs(find(goodnucleitmp),:)','color',moviecolor((imovie),:));hold on;
        plot(dcs(find(~goodnucleitmp),:)','color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('drift (\mum)');axis([0 101 0 2]);title('drift');
        subplot(2,3,4)
        plot(rmsf(:,find(goodnucleitmp)),'color',moviecolor((imovie),:));hold on;
        plot(rmsf(:,find(~goodnucleitmp)),'color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('angles');ylabel('rmsf (\mum)');axis([0 length(zrange) 0 0.4]);title('rmsf');
        subplot(2,3,5)
        bins=0:0.0025:0.15;
        rmsf1=rmsf(:,find(goodnucleitmp));
        [counts]=hist(rmsf1(:),bins);
        cumcounts=cumsum(counts)./sum(counts);
        plot(bins,cumcounts,'linewidth',2,'color',moviecolor((imovie),:));hold on;
        xlabel('rmsf(\mum)');ylabel('cummulative probability');axis([0 0.15 0 1]);title('cumrmsf');
        legend(legendids);
        subplot(2,3,6)
        plot(dr2s(:,find(goodnucleitmp)),'color',moviecolor((imovie),:));hold on;
        plot(dr2s(:,find(~goodnucleitmp)),'color',moviecolor((imovie),:),'linestyle','-.');hold on;
        xlabel('frames');ylabel('sum dr square over 6 (pixel^2)');axis([0 length(zrange) 0 1]);title('outlier');

        typermsf=[typermsf;rmsf1(:)];
    end
    name=allnames{itype};
    mkdir([verifypath,'\',name]);
    print(f1,[verifypath,'\',name,'\_existflag'],'-dpng');
    print(f2,[verifypath,'\',name,'\_zpos'],'-dpng');
    print(f3,[verifypath,'\',name,'\_drift'],'-dpng');
    print(f4,[verifypath,'\',name,'\_rmsf'],'-dpng');
    print(f5,[verifypath,'\',name,'\_cumsumrmsf'],'-dpng');
    print(f7,[verifypath,'\',name,'\_sumdr2'],'-dpng');
    print(f8,[verifypath,'\',name],'-dpng');
    
%     print(f2,[verifypath,'\zpos_',name],'-dpng');
%     print(f3,[verifypath,'\drift_',name],'-dpng');
%     print(f4,[verifypath,'\rmsf_',name],'-dpng');
%     print(f5,[verifypath,'\cumsumrmsf_',name],'-dpng');
%     print(f7,[verifypath,'\sumdr2_',name],'-dpng');
%     
    figure(f6)
    bins=0:0.0025:0.15;
    [counts]=hist(typermsf(:),bins);
    cumcounts=cumsum(counts)./sum(counts);
%     if ~isempty(strfind(name,'MBC'))
        plot(bins,cumcounts,'linestyle','-.','color',allcolors(itype,:));hold on;
%     else
%         plot(bins,cumcounts,'color',colorAll(itype,:));hold on;
%     end
    xlabel('rmsf(\mum)');ylabel('cummulative probability');axis([0 0.15 0 1])
    legend(allnames);
end
    FigureFormat(f6);
%     print(f6,[verifypath,'\allcumsumrmsf'],'-dpng');
    print(f6,[verifypath,'\','_allcumsumrmsf'],'-dpng');

    save(fullfile(rootpath,'goodnuclei.mat'),'gnmovienames','goodnuclei');
