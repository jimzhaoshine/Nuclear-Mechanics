% plot data
clear all; close all;clc;
run setup_header5.m;
rootpath='C:\nuclei\post analysis result_0.2';
load([rootpath,'\strainall.mat']);
load([rootpath,'\goodnuclei.mat']);
resultpath=[rootpath,'\result'];
mkdir(resultpath)
num_names=length(nameAll);
sizeth=1.2;
%% number of nuclei collected
figure('Position',[0 0 1200 800]);
for i=1
    meandata=zeros(num_names,1);
    sedata=zeros(num_names,1);
    stddata=zeros(num_names,1);
    for itype=1:length(nameAll2)
        meandata(itype)=length(strainall(itype).nuclei);
        B= barh(num_names+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        set(get(B,'child'),'facea',.3)
    end
    meandata1=meandata;
    for itype=1:length(nameAll2)
        if ~isempty(strainall(itype).nuclei)
            meandata(itype)=sum([strainall(itype).nuclei.good]);
            barh(num_names+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    meandata2=meandata;
    title(' Nuclei Collected (dark: good data,  light: bad data)');
    axis([0 120 0 num_names+1]);
    set(gca,'YTick',1:num_names,'YTicklabel',nameAll2(num_names:-1:1));
    xlabel('number');
    FigureFormat(gcf)
    print(gcf,[resultpath,'\number'],'-dpng');
    print(gcf,[resultpath,'\number'],'-deps');
    savefig(gcf,[resultpath,'\number']);
end
%% rmsf , cum rmsf plot for all,
figure('Position',[0 0 1200 800]);
for i=1
    mkdir([resultpath,'\rmsf']);
    mkdir([resultpath,'\cumrmsf']);
    hold off;
    zrange=find(abs(points(:,3))<0.5);
    xval=0.00:0.0025:0.1;
    rmsfcounts=zeros(length(strainall),length(xval));
    rmsfcumcounts=zeros(length(strainall),length(xval));
    for itype=1:length(strainall)
        str=strainall(itype);
        rmsftmp=[];
        for inuc=1:length(str.nuclei)
            if str.nuclei(inuc).good==1 & str.nuclei(inuc).size<sizeth
                rmsftmp=[rmsftmp,str.nuclei(inuc).rmsf(zrange)];
            end
        end
        rmsfctmp=hist(rmsftmp,xval);
        rmsfcounts(itype,:)=rmsfctmp/sum(rmsfctmp);
        rmsfcumcounts(itype,:)=cumsum(rmsfcounts(itype,:));
    end
    
    for j=1:num_names/2
        rmsfchose=[j j+num_names/2];
        plot(xval'*ones(size(rmsfchose)), rmsfcounts(rmsfchose,:)');
        legend(nameAll2(rmsfchose));
        xlabel('rmsf (\mum)');ylabel('percentage');title(catname(nameAll2(rmsfchose)));
        axis([min(xval) max(xval) 0 0.3]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\rmsf\',catname(nameAll2(rmsfchose))],'-dpng');
        print(gcf,[resultpath,'\rmsf\',catname(nameAll2(rmsfchose))],'-deps');
        savefig([resultpath,'\rmsf\',catname(nameAll2(rmsfchose))]);
    end
    
    for j=1:num_names/2
        rmsfchose=[j j+num_names/2];
        plot(xval'*ones(size(rmsfchose)), rmsfcumcounts(rmsfchose,:)');
        legend(nameAll2(rmsfchose));
        xlabel('rmsf (\mum)');ylabel('percentage');title(catname(nameAll2(rmsfchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))],'-dpng');
        print(gcf,[resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))],'-deps');
        savefig([resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))]);
    end
    
    rmsfchose=1:num_names/2;
    for ir=1:length(rmsfchose)
        plot(xval', rmsfcumcounts(rmsfchose(ir),:)','color',colorAll(rmsfchose(ir),:));hold on;
    end
    hold off;    legend(nameAll2(rmsfchose));
    xlabel('rmsf (\mum)');ylabel('percentage');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\rmsf\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\rmsf\',catname(nameAll2(rmsfchose))],'-deps');
    savefig([resultpath,'\rmsf\',catname(nameAll2(rmsfchose))]);
    
    rmsfchose=1:num_names/2;
    for ir=1:length(rmsfchose)
        plot(xval', rmsfcumcounts(rmsfchose(ir),:)','color',colorAll(rmsfchose(ir),:));hold on;
    end
    hold off;
    legend(nameAll2(rmsfchose));
    xlabel('rmsf (\mum)');ylabel('percentage');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))],'-deps');
    savefig([resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))]);
    
    rmsfchose=num_names/2+1:num_names;
    for ir=1:length(rmsfchose)
        plot(xval', rmsfcumcounts(rmsfchose(ir),:)','color',colorAll(rmsfchose(ir),:));hold on;
    end
    hold off;    legend(nameAll2(rmsfchose));
    xlabel('rmsf (\mum)');ylabel('percentage');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\rmsf\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\rmsf\',catname(nameAll2(rmsfchose))],'-deps');
    savefig([resultpath,'\rmsf\',catname(nameAll2(rmsfchose))]);
    
    rmsfchose=num_names/2+1:num_names;
    for ir=1:length(rmsfchose)
        plot(xval', rmsfcumcounts(rmsfchose(ir),:)','color',colorAll(rmsfchose(ir),:));hold on;
    end
    hold off;    legend(nameAll2(rmsfchose));
    xlabel('rmsf (\mum)');ylabel('percentage');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))],'-deps');
    savefig([resultpath,'\cumrmsf\',catname(nameAll2(rmsfchose))]);
    
end
%% rmsf vs angle plot
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='RmsfAngle';
    mkdir([resultpath,'\',analyzename]);
    hold off;
    zrange=find(abs(points(:,3))<0.5);
    nbins=40;
    dx=180/nbins;
    xval=dx:dx:180;
    rmsfangle=cell(length(xval),num_names);
    %     meanrmsfangle=zeros(length(angles),num_names);
    %     stdrmsfangle=zeros(length(angles),num_names);
    for itype=1:length(strainall)
        str=strainall(itype);
        for inuc=1:length(str.nuclei)
            if str.nuclei(inuc).good==1
                celltheta=str.nuclei(inuc).orientation;
                nuctheta=acos(points(zrange,1)*cos(celltheta)+points(zrange,2)*sin(celltheta))/pi*180;
                rmsf=str.nuclei(inuc).rmsf(zrange);
                for ith=1:length(nuctheta)
                    intangle=floor(nuctheta(ith)/dx)+1;
                    if intangle==nbins+1
                        intangle=1;
                    end
                    rmsfangle{intangle,itype}=[rmsfangle{intangle,itype},rmsf(ith)];
                end
            end
        end
    end
    meanrmsfangle=cellfun(@mean,rmsfangle);
    numrmsfangle=cellfun(@length,rmsfangle);
    stdrmsfangle=cellfun(@std,rmsfangle);
    sermsfangle=stdrmsfangle./sqrt(numrmsfangle);
    
    for j=1:num_names/2
        rmsfchose=[j j+num_names/2];
        plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));hold on;
        errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
        legend(nameAll2(rmsfchose));
        xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
        axis([min(xval) max(xval) 0.02 0.06]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
        hold off;
    end
    
    rmsfchose=1:num_names/2;
    plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));
    errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
    legend(nameAll2(rmsfchose));
    xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0.02 0.06]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
    
    rmsfchose=num_names/2+1:num_names;
    plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));
    errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
    legend(nameAll2(rmsfchose));
    xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0.02 0.06]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
end
% xlswrite('rmsfangle.xls',[['angles',nameAll];num2cell([(xval-xval(1)/2)',meanrmsfangle])],'mean');
% xlswrite('rmsfangle.xls',[['angles',nameAll];num2cell([(xval-xval(1)/2)',sermsfangle])],'se');
% xlswrite('rmsfangle.xls',[['angles',nameAll];num2cell([(xval-xval(1)/2)',stdrmsfangle])],'std');
%% fluc deducted rmsf vs angle plot
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='NonFlucRmsfAngle';
    mkdir([resultpath,'\',analyzename]);
    hold off;
    zrange=find(abs(points(:,3))<0.5);
    inzrange=abs(points(:,3))<0.5;
    nbins=40;
    dx=180/nbins;
    xval=dx:dx:180;
    rmsfangle=cell(length(xval),num_names);
    %     meanrmsfangle=zeros(length(angles),num_names);
    %     stdrmsfangle=zeros(length(angles),num_names);
    for itype=1:length(strainall)
        str=strainall(itype);
        for inuc=1:length(str.nuclei)
            if str.nuclei(inuc).good==1
                celltheta=str.nuclei(inuc).orientation;
                goodindtmp=inzrange & str.nuclei(inuc).nonfluczone;
                nuctheta=acos(points(goodindtmp,1)*cos(celltheta)+points(goodindtmp,2)*sin(celltheta))/pi*180;
                rmsf=str.nuclei(inuc).rmsf(goodindtmp);
                for ith=1:length(nuctheta)
                    intangle=floor(nuctheta(ith)/dx)+1;
                    if intangle==nbins+1
                        intangle=1;
                    end
                    rmsfangle{intangle,itype}=[rmsfangle{intangle,itype},rmsf(ith)];
                end
            end
        end
    end
    meanrmsfangle=cellfun(@mean,rmsfangle);
    numrmsfangle=cellfun(@length,rmsfangle);
    stdrmsfangle=cellfun(@std,rmsfangle);
    sermsfangle=stdrmsfangle./sqrt(numrmsfangle);
    
    for j=1:num_names/2
        rmsfchose=[j j+num_names/2];
        plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));hold on;
        errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
        legend(nameAll2(rmsfchose));
        xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
        axis([min(xval) max(xval) 0.02 0.06]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
        hold off;
    end
    
    rmsfchose=1:num_names/2;
    plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));
    errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
    legend(nameAll2(rmsfchose));
    xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0.02 0.06]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
    
    rmsfchose=num_names/2+1:num_names;
    plot(xval'*ones(size(rmsfchose)), meanrmsfangle(:,rmsfchose));
    errorbar(xval'*ones(size(rmsfchose)),meanrmsfangle(:,rmsfchose),sermsfangle(:,rmsfchose));
    legend(nameAll2(rmsfchose));
    xlabel('angles');ylabel('rmsf');title(catname(nameAll2(rmsfchose)));
    axis([min(xval) max(xval) 0.02 0.06]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',catname(nameAll2(rmsfchose))]);
end
%% number fluc collected
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);
for i=1
    meandata1=zeros(num_names,1);
    meandata2=zeros(num_names,1);
    sedata=zeros(num_names,1);
    stddata=zeros(num_names,1);
    for itype=1:length(nameAll2)
        datatmp=(strainall(itype).flucs);
        gooddata=[strainall(itype).flucs.good] &[strainall(itype).flucs.inplane];
        meandata1(itype)=sum(gooddata);
        B= barh(num_names+1-itype,meandata1(itype),'facecolor',colorAll(itype,:));hold on;
        set(get(B,'child'),'facea',.3)
    end
    for itype=1:length(nameAll2)
        if ~isempty(strainall(itype).nuclei)
            datatmp=(strainall(itype).flucs);
            gooddata=[strainall(itype).flucs.good] &[strainall(itype).flucs.inplane];
            [ s_ind, ds_ind ] = SelectMtFluc(datatmp(gooddata),fparam);
            meandata2(itype)=sum(s_ind);
            barh(num_names+1-itype,meandata2(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    title(' flucs Collected (dark: good data,  light: bad data)');
    axis([0 200 0 num_names+1]);
    set(gca,'YTick',1:num_names,'YTicklabel',nameAll2(num_names:-1:1));
    xlabel('number');
    FigureFormat(gcf)
    print(gcf,[resultpath,'\flucnumber'],'-dpng');
    print(gcf,[resultpath,'\flucnumber'],'-deps');
    savefig(gcf,[resultpath,'\flucnumber']);
end
xlswrite('number of flucutations.xls',[['value',nameAll];...
    {'total flucs'},num2cell(meandata1');... ]);
    {'MT flucs'},num2cell(meandata2');...
    ]);
%% mean std map
figure('Position',[0 0 1200 800]);
for i=[]
    mkdir([resultpath,'\MeanStdImg']);
    hold off;
    %     stdimgs=zeros(19,);
    for itype=1:length(strainall)
        str=strainall(itype);
        stdtmp=0;
        numgood=0;
        for inuc=1:length(str.nuclei)
            if str.nuclei(inuc).good==1
                stdtmp=stdtmp+str.nuclei(inuc).stdimg;
                numgood=numgood+1;
            end
        end
        stdtmp=stdtmp/numgood;
        imagesc(stdtmp,[0.02 0.05]);colorbar;colormap jet;axis image;
        title(['standard deviation map of ',str.name]);
        xlabel('longitudinal angle')
        ylabel('latitudinal angle')
        set(gca,'XTick',(0:16:64));
        set(gca, 'xticklabel', '0 | p/2 | p | 3p/2 | 2p', 'fontname', 'symbol');
        set(gca,'YTick',([2,10,18]));
        set(gca, 'yticklabel', '-p/4 | 0 | p/4', 'fontname', 'symbol');
        FigureFormat(gcf)
%         pause;
    end
    %         print(gcf,[resultpath,'\rmsfangle\',catname(nameAll2(rmsfchose))],'-dpng');
    %         print(gcf,[resultpath,'\rmsfangle\',catname(nameAll2(rmsfchose))],'-deps');
    %          savefig(gcf,[resultpath,'\rmsfangle\',catname(nameAll2(rmsfchose))]);
end
%% nuclei radius size
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='size';
    xname='radius (\mum)';
    xval=0.5:0.05:1.8;
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=[strainall(itype).nuclei.size];
            goodtmp=find([strainall(itype).nuclei.good]&[strainall(itype).nuclei.size]<sizeth);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    title(analyzename);
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    clf
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%%  cell width
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='cellwidth';
    xname='width (\mum)';
    xval=2:0.1:5;
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=[strainall(itype).nuclei.cellwidth];
            goodtmp=find([strainall(itype).nuclei.good]);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    title(analyzename);
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    clf
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%%  cell length
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='celllength';
    xname='length (\mum)';
    xval=6:0.2:12;
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=[strainall(itype).nuclei.celllength];
            goodtmp=find([strainall(itype).nuclei.good]);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    title(analyzename);
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%%  cell volume
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='cellvolume';
    xname='volume (\mum^3)';
    xval=30:2:150;
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=[strainall(itype).nuclei.cellvolume];
            goodtmp=find([strainall(itype).nuclei.good]);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    title(analyzename);
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    clf
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%% nuclei cell ratio
figure('Position',[0 0 1200 800]);
for i=1
    analyzename='ncratio';
    xname='N/C ratio (\mum^3)';
    xval=0:0.005:0.15;
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=[strainall(itype).nuclei.nucleusvolume]./[strainall(itype).nuclei.cellvolume];
            goodtmp=find([strainall(itype).nuclei.good]);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    title(analyzename);
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%% fluc location and max height angle
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0,'durationlb',0,'durationub',inf,...
    'anglelb',0,'angleub',90);
for i=1
    analyzename='maxheight_angle';
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam);
            fluctmp=fluctmp(find(s_ind));
            datatmp1=acos(cos([fluctmp.avglat]).*cos([fluctmp.avglong]))/pi*180;
            datatmp2=[fluctmp.maxheight];
            goodtmp=find([fluctmp.good]);
            datatmp1=datatmp1(goodtmp);
            datatmp2=datatmp2(goodtmp);
            xval=0:9:180;
            histdata=hist(datatmp1,xval);
            histdata=histdata/sum(histdata);
            bh=bar(xval,histdata,'r');axis([0 180 0 0.3]);hold on;
            set(get(bh,'children'),'facea',.2);
            scatter(datatmp1,datatmp2);hold off;xlabel('angle');ylabel('max height(\mum)')
            title('fluctuations distribution and max height angle plot');
            FigureFormat(gcf);
            print(gcf,[resultpath,'\',analyzename,'\',straintmp.name],'-dpng');
        end
    end
end
%% fluc location and mean height angle
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0,'durationlb',20,'durationub',200,...
    'anglelb',0,'angleub',90);
fparam2=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);
for i=1
    analyzename='meanheight_angle';
    titlestr={'fluctuations distribution and mean height angle plot','blue: non MT fluc, red: MT fluc'};
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam);
            fluctmp=fluctmp(find(s_ind));
            datatmp1=acos(cos([fluctmp.avglat]).*cos([fluctmp.avglong]))/pi*180;
            datatmp2=[fluctmp.meanheight];
            goodtmp=find([fluctmp.good]);
            datatmp1=datatmp1(goodtmp);
            datatmp2=datatmp2(goodtmp);
            xval=0:9:180;
            histdata=hist(datatmp1,xval);
            histdata=histdata/sum(histdata);
            bh=bar(xval,histdata,'b');axis([0 180 0 0.3]);hold on;
            set(get(bh,'children'),'facea',.2);
            scatter(datatmp1,datatmp2,'b');hold on;xlabel('angle');ylabel('mean height(\mum)')
            
            fluctmp=straintmp.flucs;
            [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam2);
            fluctmp=fluctmp(find(s_ind));
            datatmp1=acos(cos([fluctmp.avglat]).*cos([fluctmp.avglong]))/pi*180;
            datatmp2=[fluctmp.meanheight];
            goodtmp=find([fluctmp.good]);
            datatmp1=datatmp1(goodtmp);
            datatmp2=datatmp2(goodtmp);
            xval=0:num_names/2+1:180;
            histdata=hist(datatmp1,xval);
            histdata=histdata/sum(histdata);
            bh=bar(xval,histdata,'r');axis([0 180 0 0.3]);hold on;
            set(get(bh,'children'),'facea',.2);
            scatter(datatmp1,datatmp2,'r');hold off;xlabel('angle');ylabel('mean height(\mum)')
            
            title(titlestr);
            FigureFormat(gcf);
            print(gcf,[resultpath,'\',analyzename,'\',straintmp.name],'-dpng');
            print(gcf,[resultpath,'\',analyzename,'\',straintmp.name],'-deps');
            savefig(gcf,[resultpath,'\',analyzename,'\',straintmp.name]);
        end
    end
end

%% fluc distribution
figure('Position',[0 0 1200 800]);
% fparam=struct('maxheightlb',0,'meanheightlb',0,'durationlb',0,'durationub',inf,...
%     'anglelb',0,'angleub',90);
for i=[]
    analyzename='fluclocation';
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.flucs)
            fluctmp=straintmp.flucs;
            %             fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc( fluctmp,fparam);
            datatmp1=[fluctmp.avglong]/pi*180;
            datatmp2=[fluctmp.avglat]/pi*180;
            scatter(datatmp1(s_ind),datatmp2(s_ind),'x','linewidth',2,'markerfacecolor',colorAll(itype,:));hold on;
            scatter(datatmp1(ds_ind),datatmp2(ds_ind),'.','linewidth',2,'markerfacecolor',colorAll(itype,:));hold on;
            %             FigureFormat(gcf);
            %             print(gcf,[resultpath,'\',analyzename,'\',straintmp.name],'-dpng');
        end
    end
    xlabel('longitude');ylabel('latitude');axis([0 180 -60 60]);axis equal;
    title('fluctuations location');
end
%% number fluc per cell per min
figure('Position',[0 0 1200 800]);
fparam1=struct('maxheightlb',0,'meanheightlb',0,'durationlb',20,'durationub',200,...
    'anglelb',0,'angleub',90);
fparam2=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);
for i=1
    analyzename='flucpermin';
    xname='flucutatioons per nuclei per minute(s^{-1})';
    titlestr='total fluc number (dark and light) and MT induced fluc number (dark)';
    xval=(0:1:4)/0.5/250*60;
    mkdir([resultpath,'\',analyzename]);
    datalength=num_names;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata1=zeros(datalength,1);
    sedata1=zeros(datalength,1);
    stddata1=zeros(datalength,1);
    meandata2=zeros(datalength,1);
    sedata2=zeros(datalength,1);
    stddata2=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=zeros(1,length(straintmp.nuclei));
            for inucs=1:length(straintmp.nuclei)
                fluctmp= straintmp.nuclei(inucs).flucs;
                if ~isempty(fluctmp)
                    [ s_ind, ds_ind ] = SelectMtFluc(fluctmp([fluctmp.inplane]==1),fparam2);
                else
                    s_ind=[];
                end
                datatmp(inucs)=sum(s_ind)/0.5/250*60;
            end
            goodtmp=find([strainall(itype).nuclei.good]);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata1(itype)=mean(datatmp);
            stddata1(itype)=std(datatmp);
            sedata1(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata1(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata1,datalength:-1:1,sedata1,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            datatmp=zeros(1,length(straintmp.nuclei));
            for inucs=1:length(straintmp.nuclei)
                fluctmp= straintmp.nuclei(inucs).flucs;
                if ~isempty(fluctmp)
                    [ s_ind, ds_ind ] = SelectMtFluc( fluctmp([fluctmp.inplane]==1),fparam1);
                else
                    s_ind=[];
                end
                datatmp(inucs)=sum(s_ind)/0.5/250*60;
            end
            goodtmp=find([strainall(itype).nuclei.good]);
            datatmp=datatmp(goodtmp);
            histdata(itype,:)=hist(datatmp,xval);
            histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
            cumhistdata(itype,:)=cumsum(histdata(itype,:));
            meandata2(itype)=mean(datatmp);
            stddata2(itype)=std(datatmp);
            sedata2(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata2(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata2,datalength:-1:1,sedata2,'k');hold on;
    
    title(titlestr);
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
xlswrite('fpm.xls',[['value',nameAll];...
    {'MT mean'},num2cell(meandata1');... ]);
    {'MT std'},num2cell(stddata1');{'MT standard error'},num2cell(sedata1');...
    {'all mean'},num2cell(meandata2');...
    {'all std'},num2cell(stddata2');{'all standard error'},num2cell(sedata2');...
    ]);
%% MT rise fall
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);
for i=1
    analyzename='MTFlucTime';
    xname='time (s)';
    xval=0:10:150;
    mkdir([resultpath,'\',analyzename]);
    
    datalength=num_names;
    histdata1=zeros(datalength,length(xval));
    cumhistdata1=zeros(datalength,length(xval));
    histdata2=zeros(datalength,length(xval));
    cumhistdata2=zeros(datalength,length(xval));
    
    datalength=num_names/2;
    meandata1=zeros(datalength,1);
    sedata1=zeros(datalength,1);
    stddata1=zeros(datalength,1);
    meandata2=zeros(datalength,1);
    sedata2=zeros(datalength,1);
    stddata2=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.risetime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype,:)=hist(datatmp,xval);
            histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
            cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
            meandata1(itype)=mean(datatmp);
            stddata1(itype)=std(datatmp);
            sedata1(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata1(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata1,datalength:-1:1,sedata1,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.falltime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype,:)=hist(datatmp,xval);
            histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
            cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
            meandata2(itype)=mean(datatmp);
            stddata2(itype)=std(datatmp);
            sedata2(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata2(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata2,datalength:-1:1,sedata2,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    datalength=num_names/2;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.risetime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype+num_names/2,:)=hist(datatmp,xval);
            histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
            cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.falltime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype+num_names/2,:)=hist(datatmp,xval);
            histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
            cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), histdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
% xlswrite('mtrf.xls',[['value',nameAll(1:num_names/2)];{'rise mean'},num2cell(meandata1');...
%     {'rise std'},num2cell(stddata1');{'rise standard error'},num2cell(sedata1');...
%    {'fall mean'},num2cell(meandata2');... ]);
%     {'fall std'},num2cell(stddata2');{'fall standard error'},num2cell(sedata2');...
%     ]);
%% non MT rise fall
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',20,'durationub',200,...
    'anglelb',30,'angleub',80);
for i=1
    analyzename='NonMTFlucTime';
    xname='time (s)';
    xval=0:10:150;
    mkdir([resultpath,'\',analyzename]);
    
    datalength=num_names;
    histdata1=zeros(datalength,length(xval));
    cumhistdata1=zeros(datalength,length(xval));
    histdata2=zeros(datalength,length(xval));
    cumhistdata2=zeros(datalength,length(xval));
    
    datalength=num_names/2;
    meandata1=zeros(datalength,1);
    sedata1=zeros(datalength,1);
    stddata1=zeros(datalength,1);
    meandata2=zeros(datalength,1);
    sedata2=zeros(datalength,1);
    stddata2=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.risetime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype,:)=hist(datatmp,xval);
            histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
            cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
            meandata1(itype)=mean(datatmp);
            stddata1(itype)=std(datatmp);
            sedata1(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata1(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata1,datalength:-1:1,sedata1,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp(abs([fluctmp.avglat])<=30/180*pi);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.falltime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype,:)=hist(datatmp,xval);
            histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
            cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
            meandata2(itype)=mean(datatmp);
            stddata2(itype)=std(datatmp);
            sedata2(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata2(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata2,datalength:-1:1,sedata2,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    datalength=num_names/2;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.risetime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype+num_names/2,:)=hist(datatmp,xval);
            histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
            cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.falltime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&ds_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype+num_names/2,:)=hist(datatmp,xval);
            histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
            cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), histdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%% MT rise fall 2
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.08,'durationlb',40,'durationub',100,...
    'anglelb',30,'angleub',80);
for i=1
    analyzename='MTFlucTime2';
    xname='time (s)';
    xval=0:10:150;
    mkdir([resultpath,'\',analyzename]);
    
    datalength=num_names;
    histdata1=zeros(datalength,length(xval));
    cumhistdata1=zeros(datalength,length(xval));
    histdata2=zeros(datalength,length(xval));
    cumhistdata2=zeros(datalength,length(xval));
    
    datalength=num_names/2;
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.risetime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype,:)=hist(datatmp,xval);
            histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
            cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.falltime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype,:)=hist(datatmp,xval);
            histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
            cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    datalength=num_names/2;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.risetime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype+num_names/2,:)=hist(datatmp,xval);
            histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
            cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.falltime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype+num_names/2,:)=hist(datatmp,xval);
            histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
            cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), histdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end

%% all rise fall
figure('Position',[0 0 1200 800]);
fparam=struct('maxheightlb',0,'meanheightlb',0.0,'durationlb',20,'durationub',200,...
    'anglelb',0,'angleub',90);
for i=1
    analyzename='AllFlucTime';
    xname='time (s)';
    xval=0:10:150;
    mkdir([resultpath,'\',analyzename]);
    
    datalength=num_names;
    histdata1=zeros(datalength,length(xval));
    cumhistdata1=zeros(datalength,length(xval));
    histdata2=zeros(datalength,length(xval));
    cumhistdata2=zeros(datalength,length(xval));
    
    datalength=num_names/2;
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.risetime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'rise')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype,:)=hist(datatmp,xval);
            histdata1(itype,:)=histdata1(itype,:)/sum(histdata1(itype,:));
            cumhistdata1(itype,:)=cumsum(histdata1(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype);
        if ~isempty(straintmp.nuclei)
            fluctmp=straintmp.flucs;
            fluctmp=fluctmp([fluctmp.inplane]==1);
            [ s_ind, ds_ind ] = SelectMtFluc(fluctmp,fparam);
            datatmp=[fluctmp.risetime];
            goodtmp=find([fluctmp.good]&~strcmp({fluctmp.noisydata},'fall')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype,:)=hist(datatmp,xval);
            histdata2(itype,:)=histdata2(itype,:)/sum(histdata2(itype,:));
            cumhistdata2(itype,:)=cumsum(histdata2(itype,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename]);
    hold off;
    
    datalength=num_names/2;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.risetime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'rise')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata1(itype+num_names/2,:)=hist(datatmp,xval);
            histdata1(itype+num_names/2,:)=histdata1(itype+num_names/2,:)/sum(histdata1(itype+num_names/2,:));
            cumhistdata1(itype+num_names/2,:)=cumsum(histdata1(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
            set(get(B,'child'),'facea',.3);
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    for itype=1:datalength
        straintmp=strainall(itype+num_names/2);
        if ~isempty(straintmp.nuclei)
            [ s_ind, ds_ind ] = SelectMtFluc( straintmp.flucs,fparam);
            datatmp=[straintmp.flucs.falltime];
            goodtmp=find([straintmp.flucs.good]&~strcmp({straintmp.flucs.noisydata},'fall')&s_ind);
            datatmp=datatmp(goodtmp);
            histdata2(itype+num_names/2,:)=hist(datatmp,xval);
            histdata2(itype+num_names/2,:)=histdata2(itype+num_names/2,:)/sum(histdata2(itype+num_names/2,:));
            cumhistdata2(itype+num_names/2,:)=cumsum(histdata2(itype+num_names/2,:));
            meandata(itype)=mean(datatmp);
            stddata(itype)=std(datatmp);
            sedata(itype)=std(datatmp)./sqrt(length(datatmp));
            B=barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
        end
    end
    %     herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    
    title('MT induced fluc: rise time (light+dark) and fall time (dark)');
    axis([min(xval) max(xval) 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll2((datalength:-1:1)+num_names/2));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC'],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\',analyzename,'MBC']);
    hold off;
    
    
    clf
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), histdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), histdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    end
    
    for jj=1:num_names/2
        repchose=jj;%[jj jj+num_names/2];
        plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)','r');hold on;
        plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','b');hold off;
        legend([nameAll2{repchose},'rise'],[nameAll2{repchose},'fall']);
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 1]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
        print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
        savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    end
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), histdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), histdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\hist',catname(nameAll2(repchose))]);
    
    repchose=1:num_names/2;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
    
    repchose=num_names/2+1:num_names;
    plot(xval'*ones(size(repchose)), cumhistdata1(repchose,:)');hold on;
    plot(xval'*ones(size(repchose)), cumhistdata2(repchose,:)','-.');hold off;
    title('solid: rise   dash: fall');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-dpng');
    print(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))],'-deps');
    savefig(gcf,[resultpath,'\',analyzename,'\cumhist',catname(nameAll2(repchose))]);
end
%%
% fparam=struct('maxheightlb',0.1,'meanheightlb',0,'durationlb',20,'durationub',200,...
%     'anglelb',0,'angleub',90);
