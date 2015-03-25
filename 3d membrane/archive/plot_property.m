function plot_property(analyzename,xname,xval,fun,funfilter) 

% analyzename='size';
%     xname='radius (\mum)';
%     xval=0:1;
    datalength=16;
    histdata=zeros(datalength,length(xval));
    cumhistdata=zeros(datalength,length(xval));
    meandata=zeros(datalength,1);
    sedata=zeros(datalength,1);
    stddata=zeros(datalength,1);
    figure(1001)
    for itype=1:datalength
        straintmp=strainall(itype);
%         datatmp=[straintmp.nuclei.size];
%         goodtmp=find([straintmp.nuclei.good]);
        datatmp=fun(straintmp);
        goodtmp=funfilter(straintmp);
        datatmp=datatmp(goodtmp);
        histdata(itype,:)=hist(datatmp,xval);
        histdata(itype,:)=histdata(itype,:)/sum(histdata(itype,:));
        cumhistdata(itype,:)=cumsum(histdata(itype,:));
        meandata(itype)=mean(datatmp);
        stddata(itype)=std(datatmp);
        sedata(itype)=std(datatmp)./sqrt(length(datatmp));
        barh(datalength+1-itype,meandata(itype),'facecolor',colorAll(itype,:));hold on;
    end
    
    herrorbar(meandata,datalength:-1:1,stddata,'k');hold on;
    herrorbar(meandata,datalength:-1:1,sedata,'k');hold on;
    title(analyzename);
    axis([1 1.4 0 datalength+1]);
    set(gca,'YTick',1:datalength,'YTicklabel',nameAll(datalength:-1:1));
    xlabel(xname);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',analyzename],'-dpng');
    for jj=1:8
        repchose=[jj jj+8];
        plot(xval'*ones(size(repchose)), histdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll(repchose))],'-dpng');
    end
    
    for jj=1:8
        repchose=[jj jj+8];
        plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
        legend(nameAll2(repchose));
        xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
        axis([min(xval) max(xval) 0 .4]);
        FigureFormat(gcf)
        print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll(repchose))],'-dpng');
    end
    repchose=1:8;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll(repchose))],'-dpng');
    
    repchose=9:16;
    plot(xval'*ones(size(repchose)), histdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 .4]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll(repchose))],'-dpng');
    
    repchose=1:8;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\',analyzename,'\',catname(nameAll(repchose))],'-dpng');
    
    repchose=9:16;
    plot(xval'*ones(size(repchose)), cumhistdata(repchose,:)');
    legend(nameAll2(repchose));
    xlabel(xname);ylabel('percentage');title(catname(nameAll2(repchose)));
    axis([min(xval) max(xval) 0 1]);
    FigureFormat(gcf)
    print(gcf,[resultpath,'\,',analyzename,'\',catname(nameAll(repchose))],'-dpng');
end