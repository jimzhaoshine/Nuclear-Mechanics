%find cell orientation
close all;
clear all;
rootpath1='C:\nuclei';
datapath1=[rootpath1,'\data'];
folders1=dir(datapath1);
numfolders=length(folders1);

%% find orientation by hand
for ifolder=12:numfolders
    strpath=folders1(ifolder).name;
    if regexpi(strpath,'wild type|heh1|ima1|heh2')
        display(['Process folder ',strpath]);
        moviename=dir([datapath1,'\',strpath,'\*.mat']);
        if ~isempty(moviename)
            %%
            for jmovie=1:length(moviename)
                %%
                load([datapath1,'\',strpath,'\',moviename(jmovie).name]);
                nm.path=[datapath1,'\',strpath];
%                 nm.filein=strrep([datapath,'\',strpath,'\',moviename(j).name],'.mat','.dv');
                nm.loadmovie;
                display(['Process movie ',moviename(jmovie).name]);
                %%
                iframe=1;
                img=nm.proj(iframe);
                wsize=50;
                
                x=zeros(nm.num_nuc,4);
                y=zeros(nm.num_nuc,4);

%                 orientation=zeros(nm.num_nuc,1);
%                 celllength=zeros(nm.num_nuc,1);
%                 cellwidth=zeros(nm.num_nuc,1);

                for inuc=1:size(nm.nuclei,2)
                    nuc=nm.nuclei{iframe,inuc};
                    cnt=round(nuc.center);
                    wimg0=WindowImageUS(cnt(1),cnt(2),wsize,img);
                    clf
                    SI(wimg0);hold on;
                    for ichoose=1:4
                        if ichoose<=2
                            title('choose long axis');
                        else
                            title('choose short axis')
                        end
                        [x(inuc,ichoose),y(inuc,ichoose)]=ginput(1);
                        plot(x(inuc,ichoose),y(inuc,ichoose),'ro');
                        text(x(inuc,ichoose)+2,y(inuc,ichoose)+2,num2str(ichoose),...
                            'color','r');
                    end                    
                end
                orientation=(atan((y(:,2)-y(:,1))./(x(:,2)-x(:,1))))*180/pi;
                celllength=sqrt((y(:,2)-y(:,1)).^2+(x(:,2)-x(:,1)).^2);
                cellwidth=sqrt((y(:,4)-y(:,3)).^2+(x(:,4)-x(:,3)).^2);

                nm.orientation=orientation;
                nm.celllength=celllength;
                nm.cellwidth=cellwidth;
                nm.cpts4.x=x;
                nm.cpts4.y=y;
                %
                nm.save_contour(0);
                close all;
            end
        end
    end
end

%% save all length and all width
saveind=1;
clear allorientations;
clear allcelllengths;
clear allcellwidths;
clear allcellvols;
clear movienames1;
for ifolder=3:numfolders
    strpath=folders1(ifolder).name;
    if regexpi(strpath,'wild type|heh1|ima1|heh2')
        display(['Process folder ',strpath]);
        moviename=dir([datapath1,'\',strpath,'\*.mat']);
        if ~isempty(moviename)
            for jmovie=1:length(moviename)
                load([datapath1,'\',strpath,'\',moviename(jmovie).name]);
                
                allorientations{saveind}=nm.orientation;
                allcelllengths{saveind}=nm.celllength;
                allcellwidths{saveind}=nm.cellwidth;
                movienames1{saveind}=moviename(jmovie).name;
                saveind=saveind+1;
            end
        end
    end
end
save([datapath1,'\loc.mat'],'allorientations','allcelllengths','allcellwidths','movienames1');

