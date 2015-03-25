%% remake everything into movies
clear all
close all
rootpath='C:\nuclei';
folders=dir([rootpath,'\data']);
categoryID=0;
moviestats=[];

for ifolders=1:length(folders)
    strpath=folders(ifolders).name;
    if regexpi(strpath,'wild type|heh1|ima1|heh2')
        display(['Process folder ',strpath]);
        moviename=dir([rootpath,'\data\',strpath,'\*.mat']);
        if ~isempty(moviename)
            categoryID=categoryID+1;
            category=moviename(1).name;
            category=category(1:end-7);
%             num_total=0;
%             num_completed=0;
            for j=1:length(moviename)
                %%
                load([rootpath,'\data\',strpath,'\',moviename(j).name]);
                display([moviename(j).name,' loaded']);
%                 num_total=num_total+nm.num_nuc;
%                 if nm.continuefrom_frame>=102%nm.endframe
%                     num_completed=num_completed+nm.num_nuc;
%                     nm.path=['C',nm.path(2:9),'\data',nm.path(10:end)];
                    nm.loadmovie(10);
                    %% save nuclei selection to the folder
                    f1=figure(61114001);
                    SI(nm.grab(1,5));hold on;
                    title('nuclei selection at frame=1 z=5');axis off;
                    if ~isempty(nm.nuclei)
                    cnt=[];
                    for inuc=1:nm.num_nuc
                        cnt(inuc,:)=nm.nuclei{1,inuc}.origin_new;
                    end
                    if ~isempty(cnt)
                        plot(cnt(:,1),cnt(:,2),'o','Linewidth',2)
                    end
                    end
                    mkdir([rootpath,'\nuclei_selection']);
                    print(f1,[rootpath,'\nuclei_selection\',nm.filename],'-dpng');
                    close(f1);
                    
                    %% make movie (not drift corrected)
                    
                    if 0
                    num_submovie=ceil(nm.num_nuc/8);
                    wsize=25;
                    points=nm.points;
                    faces=nm.faces;
                    sfw=200;
                    for isub=1:num_submovie
                        fhandle=figure('Position',[0 50 1600 820],'menubar','none','color','k');
                        set(fhandle,'Renderer','zbuffer')
                        opengl('software');
                        daObj=VideoWriter([rootpath,'\analyzed_movies\',nm.filename,'_',num2str(isub)],'Uncompressed AVI');
                        %     daObj=VideoWriter([rootpath,'\analyzed_movies\',nm.filename,'_',num2str(isub)],'MPEG-4');
                        %     daObj.set('Quality',100);
                        daObj.FrameRate=1;
                        open(daObj);
                        inuc_subs=(1+(isub-1)*8):min([isub*8,nm.num_nuc]);
                        for iframe=10:nm.endframe
                            set(fhandle,'name',['frame: ',num2str(iframe)]);
                            clf
                            img=nm.grab3(iframe);
                            for isub_nuc=1:length(inuc_subs)
                                inuc=inuc_subs(isub_nuc);
                                nuc=nm.nuclei{iframe,inuc};
                                cnt=round(nuc.origin);
                                istack=cnt(3)-1:cnt(3)+1;
                                for i=1:3
                                    dx=nuc.origin(1)-round(nuc.origin(1));
                                    dy=nuc.origin(2)-round(nuc.origin(2));
                                    if istack(i)<=10 && istack(i)>0
                                        pos=[nuc.contour(istack(i)).x+dx,nuc.contour(istack(i)).y+dy]+(wsize+1);
                                        ahandle=axes('Units','Pixels','Position',[(isub_nuc-1)*sfw i*sfw-sfw sfw sfw]);
                                        wimg=WindowImageUS(cnt(1),cnt(2),wsize,squeeze(img(:,:,istack(i))));
                                        imagesc(wimg);axis image;colormap gray;hold on;
                                        plot(pos(:,1),pos(:,2),'r-','linewidth',1);
                                        text(5,5,['stack: ',num2str(istack(i))],'color','m');
                                        text(5,10,['area: ',num2str(nuc.contour(istack(i)).area*nm.pix^2,'%.2f'),...
                                            'nm^2'],'color','c');
                                        set(ahandle,'box','on','XColor',[1 1 1],'YColor',[1 1 1])
                                        set(ahandle,'xtick',[],'ytick',[],'linewidth',1)
                                    end
                                end
                                ahandle=axes('Units','Pixels','Position',[(isub_nuc-1)*sfw 4*sfw-sfw sfw sfw]);
                                pts=[nuc.r.*points(:,1),nuc.r.*points(:,2),nuc.r.*points(:,3)];
                                TR = triangulation(faces,pts);
                                trisurf(TR,'FaceColor','red','EdgeColor','none');
                                text(-15,15,25,['nuclei index: ',num2str(inuc)],'color','g','linewidth',5);
                                if isub_nuc==1
                                    text(-15,15,20,['frame: ',num2str(iframe)],'color','y','linewidth',5);
                                end
                                text(-15,15,-20,['center: ',num2str(nuc.origin(1)*nm.pix,'%.2f '),' / ',...
                                    num2str(nuc.origin(2)*nm.pix,'%.2f '),' / ',num2str(nuc.origin(3)*nm.vox,'%.2f ')],...
                                    'color','b','linewidth',5);
                                
                                axis([-15 15 -15 15 -15 15]);
                                axis off
                                view(3);
                                daspect([1 1 1/0.85])
                                camlight
                                lighting gouraud
                                set(ahandle,'box','on','XColor',[1 1 1],'YColor',[1 1 1])
                                set(ahandle,'xtick',[],'ytick',[],'linewidth',1)
                            end
                            drawnow;
                            pause(0.02);
                            frame=getframe(fhandle);
%                             writeVideo(daObj,frame);
                        end
                        close(fhandle);
                        close(daObj);
                    end
%                 else
%                     display([moviename(j).name,' is strange']);
%                 end
                    end
            end
        end
    end
end

%%
%%

