% make movies for all samples

dirpath='C:\nuclei\data';
dirpathsave='C:\nuclei\post analysis result_0.2\movie';
mkdir(dirpathsave);
names={'wild type\wild_type_04','heh1\heh1_07','heh1heh2\heh1heh2_06',...
    'heh1heh2ima1\heh1heh2ima1_01','heh1ima1\heh1ima1_04','heh2\heh2_02',...
    'heh2ima1\heh2ima1_02','ima1\ima1_03'};
index=[3,6,7,3,3,5,3,11];
%%
for i1=1:length(names)
    filename=fullfile(dirpath,[names{i1},'.mat']);
    load(filename)
    numframes=101;
    nm.loadmovie(numframes*10);
    %%
    daObj=VideoWriter(fullfile(dirpathsave,fileparts(names{i1})),'MPEG-4');
%     daObj.set('Quality',100);
    daObj.FrameRate=3;
    open(daObj);
    %%
    axesize=720;
    figure('Position',[0 50 axesize*4/3 axesize]);
    inuc=index(i1);
    wsize=25;
    nuc0=nm.nuclei{1,inuc};
    cnt0=floor(nm.nuclei{1,inuc}.origin_new);
    dcnt0=nuc0.origin_new-cnt0;
    flipstate=dcnt0>=0.5;
    z=abs(nm.points(:,3));
    capregion=z>.5;
    allr=zeros(length(nuc0.r),numframes);
    for iframe=1:numframes
        allr(:,iframe)=nm.nuclei{iframe,inuc}.r;
    end
    meanr=mean(allr,2);
    dr=allr-meanr*ones(1,numframes);
    filterR=@(dr)dr(capregion,:).*(((1-z(capregion))/.5+0.15)*ones(1,size(dr,2))).^1.3;
    %     dr(capregion)=dr(capregion).*((1-z(capregion))/.5).^2;
    dr(capregion,:)=filterR(dr);
    ranger=2*std(dr(:));
    for iframe=6:numframes-10
        clf
        nuc=nm.nuclei{iframe,inuc};
        cnt=floor(nuc.origin_new);
        dcnt=nuc.origin_new-cnt;
        % here is an assumption that change in cnt is gradual, so there is
        % time to flip before going to next integer, also requires that
        % margin is <0.4
        margin=.2;
        for idimen=1:3
            if flipstate(idimen)==0 && dcnt(idimen)>=0.5+margin
                flipstate(idimen)=1;
            elseif flipstate(idimen)==1 && dcnt(idimen)<0.5-margin
                flipstate(idimen)=0;
            end
        end
        cnt=floor(nuc.origin_new)+flipstate;
        dcnt=nuc.origin_new-cnt;
        
        dri=nuc.r-meanr;
        dri(capregion,:)=filterR(dri);
%         colorind=(dri+ranger)/ranger/2;
%         colorind=dri/ranger;
        colorind=dri*0.16/0.2;
        colorind(colorind>1)=1;
        colorind(colorind<0)=0;
        cm=jet(256);
        colorind=floor(colorind*255)+1;
        rgb=cm(colorind,:);
        points=nm.points;
        faces=nm.faces;
        pts=[nuc.r.*points(:,1),-nuc.r.*points(:,2),nuc.r.*points(:,3)];
        img=nm.grab3(iframe);
        istack=cnt(3)-1:cnt(3)+1;
        cutsize=30;
        a0=axes('Unit','Pixels','Position',[cutsize cutsize*3 axesize-2*cutsize axesize-2*cutsize]);
        p.faces=faces;
        p.vertices=pts;
%         TR = triangulation(faces,pts);
%         trisurf(TR,'FaceColor','red','EdgeColor','none');
        patch(p,'EdgeColor','none','FaceVertexCData',rgb,'FaceColor','interp');
%         colormap jet;
        axis([-15 15 -15 15 -15 15]);
        view([0 0 1]);
        daspect([1 1 1/0.85])
%         camlight
        lighting gouraud
%         box off;
        tickslabel=cellfun(@(x)[num2str(x),''],num2cell(linspace(0,200,5)),'UniformOutput',0);
        c=colorbar('Location','south','Ticks',0:0.25:1,'TickLabels',tickslabel,...
            'FontSize',20,'FontWeight','bold'); 
        colormap(a0,jet);
        text(-1.2,-15,0,'(nm)','FontSize',20,'FontWeight','bold');
        axis off;
       
        for i=1:3
            dx=nuc.origin_new(1)-cnt(1);
            dy=nuc.origin_new(2)-cnt(2);
            if istack(i)<=10 && istack(i)>0
                axes('Unit','Pixels','Position',[axesize (i-1)*axesize/3 axesize/3 axesize/3]);
                pos=[nuc.contour(istack(i)).x+dx,nuc.contour(istack(i)).y+dy]+(wsize+1);
                wimg=WindowImageUS(cnt(1),cnt(2),wsize,squeeze(img(:,:,istack(i))));
                imagesc(wimg);axis image;colormap gray;hold on;
                plot(pos(:,1),pos(:,2),'r-','linewidth',2);   
                axis off;
                rectangle('Position',[40 45 1/0.16 1],'EdgeColor','none','FaceColor','w');
%                 text(40,45-tsize,'0','color','w');
%                 text(41,43,'1\mum','color','w');
                text(3,3,['z stack ',num2str(istack(i))],'color','y','FontSize',20,'FontWeight','bold');
%                         surf(wimg'/max(wimg(:)));view([1 0 5]);axis off;
            end
        end
        writeVideo(daObj,getframe(gcf));
    end
    close(daObj);

end
