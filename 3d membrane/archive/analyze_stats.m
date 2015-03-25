%statistics of analyzed movie

rootpath='C:\nuclei';
folders=dir(rootpath);
categoryID=0;
for i=1:length(folders)
    strpath=folders(i).name;
    moviename=dir([rootpath,'\',strpath,'\*.mat']);
    if ~isempty(moviename)
        categoryID=categoryID+1;
        category=moviename(1).name;
        category=category(1:end-7);
        num_total=0;
        num_completed=0;
        for j=1:length(moviename)
            load([rootpath,'\',strpath,'\',moviename(j).name]);
            num_total=num_total+nm.num_nuc;
            if nm.continuefrom_frame>=102%nm.endframe
                num_completed=num_completed+nm.num_nuc;
            else
                moviename(j).name
            end
        end
        report(categoryID).name=category;
        report(categoryID).num_total=num_total;
        report(categoryID).num_completed=num_completed;
    end
end