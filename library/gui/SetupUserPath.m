function [  ] = SetupUserPath(  )
%set up global path for each user

computername=getenv('COMPUTERNAME');
username=getenv('USERNAME');


savedcomputers=[];
savedusers=[];
saveddatapath=[];
savedcodepath=[];
% saved names
% home pc
newcomputer='ZSHINEY-PC';
newuser='Zshiney';
newdatapath='D:\googledrive';
newcodepath='D:\dropbox\Dropbox\matlab';
savedcomputers = [savedcomputers,{newcomputer}];
savedusers = [savedusers,{newuser}];
saveddatapath = [saveddatapath,{newdatapath}];
savedcodepath = [savedcodepath,{newcodepath}];
% alienware


% sted


% office


% tweezer


for i=1:length(savedcomputers)
    if strcmp(savedcomputers{i},computername) && strcmp(savedusers{i},username)
        global datapath;
        global codepath;
        datapath = saveddatapath{i};
        codepath = savedcodepath{i};
        display(['datapath: ',datapath]);
        display(['codepath: ',codepath]);
    else
        warning('This computer is not set up yet');
    end
end

end

