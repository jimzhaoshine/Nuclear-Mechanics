function [ s_ind, ds_ind ] = SelectMtFluc(flucs,fparam)
%select microtubule induced fluctuation from given data
% input: fluctuation struct array
%        filter parameters array
% output: selected index of that array

if ~isempty(flucs)
    maxheights=[flucs.maxheight];
    meanheights=[flucs.meanheight];
    maxheightsrelative=[flucs.maxheight]./[flucs.size];
%     thetas=abs((mod([flucs.avglong]+pi,2*pi)-pi));
    theta=acos(cos([flucs.avglat]).*cos([flucs.avglong]))/pi*180;
    
    riseonly=find(strcmp({flucs.noisydata},'fall'));
    fallonly=find(strcmp({flucs.noisydata},'rise'));
    bothrf=find(strcmp({flucs.noisydata},'none'));
    duration=zeros(size([flucs.risetime]));
    duration(bothrf)=[flucs(bothrf).risetime]+[flucs(bothrf).falltime];
    duration(riseonly)=[flucs(riseonly).risetime]*2;
    duration(fallonly)=[flucs(fallonly).falltime]*2;
    
    filtered=maxheights>=fparam.maxheightlb & duration>= fparam.durationlb ...
        & duration<=fparam.durationub & meanheights>=fparam.meanheightlb...
        & abs(theta-90)>=fparam.anglelb &abs(theta-90)<=fparam.angleub;
    s_ind=(filtered);
    ds_ind=(~filtered);
else
    s_ind=[];
    ds_ind=[];
end

end

