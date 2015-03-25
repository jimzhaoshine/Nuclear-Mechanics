%statistics of analyzed movie
rootpath='C:\nuclei';
folders=dir(rootpath);
figure
for i=1:length(folders)
    strpath=folders(i).name;
%     if sum(strcmp(strpath,{'wild type' 'heh1' 'heh2' 'ima1'}))
    if sum(strfind(strpath,'MBC'))
        moviename=dir([rootpath,'\',strpath,'\*.mat']);
        if ~isempty(moviename)
            names={moviename.name};
            for j=1:length(names)
                [~,namei]=fileparts(names{j});
                filein=[rootpath,'\GoodNucleiIndex\',namei,'.txt'];
                if exist(filein,'file')
                    display(namei)
                    M = dlmread(filein);
                    cntsa=[M(:,1),512-M(:,2),5.5+zeros(size(M,1),1)];
                    load([rootpath,'\',strpath,'\',namei]);
                    cnt=nm.cnt_tmp;
                    chooseind=[];
                    cntf=[];
                    cntfl=[];
                    for k=1:size(cnt,1)
                        findk=0;
                        for l=1:size(cntsa,1)
                            if (cnt(k,1)-cntsa(l,1))^2+(cnt(k,2)-cntsa(l,2))^2 < 20^2
                                findk=1;
                            else
                            end
                        end
                        if findk==1
                            chooseind=[chooseind;k];
                            cntf=[cntf;cnt(k,:)];
                        else
                            cntfl=[cntfl;cnt(k,:)];
                        end
                    end
                    
                    nm.nuclei=nm.nuclei(:,chooseind);
                    nm.num_nuc=size(nm.nuclei,2);
                    savefile=fullfile([rootpath,'\result'],[nm.filename,'.mat']);
                    save(savefile);
                    
                    clf
                    plot(cntsa(:,1),cntsa(:,2),'b.');
                    hold on;
                    if ~isempty(cntf)
                        plot(cntf(:,1),cntf(:,2),'go');
                    end
                    hold on;
                    if ~isempty(cntfl)
                        plot(cntfl(:,1),cntfl(:,2),'ro');
                    end
                    legend('sarah','nuclei good','nuclei bad')
                    axis([0 512 0 512])
                end
            end
        end
%     else
%         moviename=dir([rootpath,'\',strpath,'\*.mat']);
%         if ~isempty(moviename)
%             names={moviename.name};
%             for j=1:length(names)
%                 [~,namei]=fileparts(names{j});
%                 filein=[rootpath,'\GoodNucleiIndex\',namei,'.txt'];
%                 if exist(filein,'file')
%                     display(namei)
%                     M = dlmread(filein);
%                     cntsa=[M(:,1),512-M(:,2),5.5+zeros(size(M,1),1)];
%                     load([rootpath,'\',strpath,'\',namei]);
%                     cnt=nm.cnt_tmp;

    end
end

