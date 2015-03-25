% function detect_fluc()
% rootpath='C:\nuclei';
% allpath=[];
% allfolder=dir([rootpath,'\data']);
% for i=3:length(allfolder)
%     folderfile=dir([rootpath,'\data\',allfolder(i).name,'\*.mat']);
%     folderpath=cellfun(@(x)fullfile([rootpath,'\data\',allfolder(i).name],x),{folderfile.name},'UniformOutput',0);
%     allpath=[allpath,folderpath];
% end
% allpath

% folderfile=dir([rootpath,'\data\ima1 MBC','\*.mat']);
% allpath=cellfun(@(x)fullfile([rootpath,'\data\ima1 MBC'],x),{folderfile.name},'UniformOutput',0);
   
rootpath='C:\nuclei\post analysis result_0.2\data';
folderfile=dir([rootpath,'\*.mat']);
allpath=cellfun(@(x)fullfile(rootpath,x),{folderfile.name},'UniformOutput',0);


for ifile=1:length(allpath)
%     singlerun(allpath{ifile})
ipath=allpath{ifile};
if ~isempty(strfind(ipath,'sp10')) || ~isempty(strfind(ipath,'wild_type'))
% if ~isempty(strfind(ipath,'wild_type'))
        data=load(ipath);
    nm=data.nm;
    nm.detect_fluctuation;
    save(ipath,'nm');

end
end
% end
% function singlerun(ipath)
%     data=load(ipath);
%     nm=data.nm;
%     nm.detect_fluctuation;
%     save(ipath,'nm');
% end