%post process all movies, correct for drift and interpolate 2d img
t_all=tic;
rootpath='C:\nuclei';
datapath=[rootpath,'\data'];
folders=dir(datapath);
categoryID=0;
moviestats=[];
% drift_param=0.05;
drift_param=0.2;

for i=1:length(folders)
    strpath=folders(i).name;
    if regexpi(strpath,'wild type MBC')%|heh1|ima1|heh2')
        display(['Process folder ',strpath]);
        moviename=dir([datapath,'\',strpath,'\*.mat']);
        if ~isempty(moviename)
            categoryID=categoryID+1;
            category=moviename(1).name;
            category=category(1:end-7);
            num_total=0;
            num_completed=0;
            for j=1:length(moviename)
                tic
                load([datapath,'\',strpath,'\',moviename(j).name]);
                display([moviename(j).name,' loaded']);
                num_total=num_total+nm.num_nuc;
                if nm.continuefrom_frame>=102%nm.endframe
                    num_completed=num_completed+nm.num_nuc;
                    display(['correcting drift for: ',moviename(j).name])
                    nm.correct_drift(.2);
                    display(['interpolating centralband for: ',moviename(j).name])
                    nm.centralband_all;
                    savefile=fullfile(['C:\nuclei\post analysis result_',num2str(drift_param),'\data'],...
                        [nm.filename,'.mat']);
                    save(savefile,'nm');
                    display(['data saved for: ',moviename(j).name])
                    toc
                else
                    display([moviename(j).name,' is strange']);
                end
            end
        end
    end
end
exp_all=toc(t_all)