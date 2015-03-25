% function process()
allpath=[];
rootpath='C:\nuclei\data';
allfolder=dir(fullfile(rootpath,'*cdc25'));
for i=1:length(allfolder)
    folderfile=dir([rootpath,'\',allfolder(i).name,'\*.mat']);
    folderpath=cellfun(@(x)fullfile(rootpath,allfolder(i).name,x),{folderfile.name},'UniformOutput',0);
    allpath=[allpath,folderpath];
end
allpath
%%
% allfile=dir('*.mat');
% allpath={allfile.name};
parfor i=1:length(allpath)
    allpath{i}
    singleworker(allpath{i});
    
end
% end



