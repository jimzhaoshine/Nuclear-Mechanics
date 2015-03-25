%% find mean fluc location relative to cell orientation
%set header
load([path,'\pt_filter\alltracks.mat']);
nameNMBC={'wild_type','heh1','heh2','ima1','heh1heh2','heh1ima1','heh2ima1','heh1heh2ima1'};
nameMBC={'sp10_MBC','heh1_MBC','heh2_MBC','ima1_MBC','heh1heh2_MBC','heh1ima1_MBC','heh2ima1_MBC','heh1heh2ima1_MBC'};
nameAll=[nameNMBC,nameMBC];

%parse index
nuctype=cell(size(nucind));
nucmovie=cell(size(nucind));
nucnumind=zeros(size(nucind));
for inuci=1:length(nucind)
    namefull=nucind{inuci};
    strtmpind=regexp(namefull,'_');
    nuctype{inuci}=namefull(1:strtmpind(end-1)-1);
    nucmovie{inuci}=namefull(1:strtmpind(end)-1);
    nucnumind(inuci)=str2double(namefull(strtmpind(end)+1:end));
end

% extract locations
allpsi=cell(size(nameAll));
alltheta=cell(size(nameAll));
for iname=1:length(nameAll)
    name1=nameAll{iname};
    typeind=cellfun(@(x)(strcmp(x,name1)),nuctype);
    typetrj=trj(typeind);
    typemovie=nucmovie(typeind);
    typenucnumind=nucnumind(typeind);
    uniquemovie=unique(typemovie);
    psitmp=[];
    thetatmp=[];
    for imovie=1:length(uniquemovie)
        display(['processing ',uniquemovie{imovie}] );
%         load(fullfile(datapath,[uniquemovie{imovie},'.mat']));
        movieind=cellfun(@(x)(strcmp(x,uniquemovie{imovie})),typemovie);
        movietrj=typetrj(movieind);
        movienucnumind=typenucnumind(movieind);
        orientmp=allorientations{cellfun(@(x)strcmp(x(1:end-4),uniquemovie{imovie}),movienames)};
        orientations=orientmp(movienucnumind);
        for inuctrj=1:length(movietrj)
            nuctrj=movietrj{inuctrj};
            for itrj=1:nuctrj(end,end)
                subtrj=nuctrj(nuctrj(:,end)==itrj,:);
                if ~isempty(subtrj)
                    if max(subtrj(:,4))>0.1
                        subtrjnonext=subtrj(subtrj(:,5)~=0,:);
                        mx=mean(subtrjnonext(:,1));
                        my=mean(subtrjnonext(:,2));
                        psi= (mx-8)/64*360-orientations(inuctrj);
                        theta= (my-8-10)/64*360;
                        psirad = psi*pi/180;
                        thetarad = theta*pi/180;
                        psi= acos(abs(cos(psirad))*cos(thetarad))/pi*180;
                        psitmp=[psitmp,psi];
                        thetatmp=[thetatmp,theta];
                    end
                end
            end
        end
    end
    allpsi{iname}=psitmp;
    alltheta{iname}=thetatmp;
end
save([path,'\pt_filter\loc_orientation.mat'],'allpsi','alltheta');

  


