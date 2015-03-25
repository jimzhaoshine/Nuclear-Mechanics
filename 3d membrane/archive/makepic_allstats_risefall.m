
%% filtering all tracks based 
clear all
path='C:\nuclei\post analysis result_0.1';
load([path,'\pt_filter\alltracks.mat']);
mhth=0.1;% filter out all the max height that is smaller than mhth
trj_new=cell(size(trj));
for i =1:length(trj)
    tmptrj1=trj{i};
    tmptrj_new=[];
        for itrj=1:tmptrj1(end,end)
            tmptrj2=tmptrj1(tmptrj1(:,end)==itrj,1:end);
            if ~isempty(tmptrj2)
                if max(tmptrj2(:,4))>mhth
                    tmptrj_new=[tmptrj_new;tmptrj2];
                end
            end
        end
        trj_new{i}=tmptrj_new;
end

nucind=nucind(cellfun(@(x)~isempty(x),trj_new));
trj=trj_new(cellfun(@(x)~isempty(x),trj_new));


savepicdir=[path,'\pic_',num2str(mhth)];
if ~exist(savepicdir,'dir')
    mkdir(savepicdir);
end
%% group all trajectories together

truncated_nucind=cell(size(nucind));
for ii=1:length(nucind)
    namefull=nucind{ii};
    strtmpind=regexp(namefull,'_');
    name=namefull(1:strtmpind(end-1)-1);
    truncated_nucind{ii}=name;
end
mutant_cat=unique(truncated_nucind);
% save all the tracjetories together
alltracks0=cell(size(mutant_cat));
for icat=1:length(mutant_cat)
    subtrj=trj(strcmp(truncated_nucind,mutant_cat{icat}));
    subtracks=[];
    subtrackid=1;
    for isub=1:length(subtrj)
        tmptrj1=subtrj{isub};
        for itrj=1:tmptrj1(end,end)
            tmptrj2=tmptrj1(tmptrj1(:,end)==itrj,1:end);
            if ~isempty(tmptrj2)
                subtracks{subtrackid}=tmptrj2;
                subtrackid=subtrackid+1;
            end
        end
    end
    alltracks0{icat}=subtracks;
end
nameNMBC={'wild_type','heh1','heh2','ima1','heh1heh2','heh1ima1','heh2ima1','heh1heh2ima1'};
orderNMBC=cellfun(@(x)find(strcmp(mutant_cat,x)),nameNMBC);
nameMBC={'sp10_MBC','heh1_MBC','heh2_MBC','ima1_MBC','heh1heh2_MBC','heh1ima1_MBC','heh2ima1_MBC','heh1heh2ima1_MBC'};
orderMBC=cellfun(@(x)find(strcmp(mutant_cat,x)),nameMBC);
nameIntercleave=reshape([nameNMBC;nameMBC],1,16);
nameAll=[nameNMBC,nameMBC];
numM=length(nameNMBC);
alltracks=alltracks0([orderNMBC,orderMBC]);
%% single value comaprision across all mutant with errorbar
%stats on how many nuclei and how many fluction detected for each mutant
numNuc=zeros(2*numM,1);
%calculate or load num of nuclei
if ~exist([path,'\pt_filter\stat.mat'],'file')
    files=dir([path,'\data']);
    for ifile=3:length(files);
        filename=files(ifile).name;
        [~,name,ext]=fileparts(filename);
        display(['processing ',name]);
        if strcmp(ext,'.mat')
            load(fullfile([path,'\data'],filename));
            numNuc(strcmp(name(1:end-3),nameAll))=...
                numNuc(strcmp(name(1:end-3),nameAll))+size(nm.nuclei,2);
        end
    end
    save([path,'\pt_filter\stat.mat'],'numNuc');
else load([path,'\pt_filter\stat.mat']);
end
%% stats on number, and totaly length of fluctuation, and max peak size, on each cell
numFluc_pn=cell(1,numM*2);
lengthFluc_pn=cell(1,numM*2);
maxheights=cell(1,numM*2);
durations=cell(1,numM*2);
meanwidth=cell(1,numM*2);
for icat=1:numM*2
    subtrj=trj(strcmp(truncated_nucind,nameAll{icat}));
    subnumFluc=zeros(size(subtrj));
    sublengthFluc=zeros(size(subtrj));
    submaxheight=[];
    subduration=[];
    submeanwidth=[];
    for isub=1:length(subtrj)
        tmptrj1=subtrj{isub};
        subflucnum=0;
        subfluclength=0;
        for itrj=1:tmptrj1(end,end)
            tmptrj2=tmptrj1(tmptrj1(:,end)==itrj,1:end);
            if ~isempty(tmptrj2)
                subflucnum=subflucnum+1;
                subfluclength=subfluclength+sum(tmptrj2(:,5)~=0);
                submaxheight=[submaxheight,max(tmptrj2(:,4))];
                subduration=[subduration,sum(tmptrj2(:,5)~=0)];
                submeanwidth=[submeanwidth,mean(tmptrj2(tmptrj2(:,5)~=0,5))];
            end
        end
        subnumFluc(isub)=subflucnum;
        sublengthFluc(isub)=subfluclength;
    end
    numFluc_pn{icat}=[subnumFluc;zeros(numNuc(icat)-length(subnumFluc),1)];
    lengthFluc_pn{icat}=[sublengthFluc;zeros(numNuc(icat)-length(sublengthFluc),1)];
    maxheights{icat}=submaxheight;
    durations{icat}=subduration*2.5;
    meanwidth{icat}=submeanwidth/64*360;
end

% average and error on average
avg_numFluc=cellfun(@mean,numFluc_pn);
error_numFluc=cellfun(@std,numFluc_pn)./sqrt(numNuc)';
avg_lengthFluc=cellfun(@mean,lengthFluc_pn);
error_lengthFluc=cellfun(@std,lengthFluc_pn)./sqrt(numNuc)';
avg_maxheights=cellfun(@mean,maxheights);
error_maxheights=cellfun(@std,maxheights);%./sqrt(cellfun(@length,maxheights));
avg_durations=cellfun(@mean,durations);
error_durations=cellfun(@std,durations);%./sqrt(cellfun(@length,durations));
avg_meanwidth=cellfun(@mean,meanwidth);
error_meanwidth=cellfun(@std,meanwidth);%./sqrt(cellfun(@length,meanwidth));

% average number of fluctuation
figure(1001)
barh(18-2*(1:numM),avg_numFluc(1:8),0.5,'r'); hold on;
barh(17-2*(1:numM),avg_numFluc(9:16),0.5,'b');
herrorbar(avg_numFluc,[(16:-2:2),15:-2:1],error_numFluc,'k');
set(gca,'YTick',1:16,'YTicklabel',nameIntercleave(end:-1:1));
title('average number of fluctuations detected per nuclei (errorbar standard_error)')
xlabel('average number')
print(gcf,[savepicdir,'\fluc_per_nuc'],'-dpng');
% print(gcf,[savepicdir,'\fluc_per_nuc'],'-deps');

% average sumed length of fluctuations
figure(1002)
barh(18-2*(1:numM),avg_lengthFluc(1:8),0.5,'r'); hold on;
barh(17-2*(1:numM),avg_lengthFluc(9:16),0.5,'b');hold on;
herrorbar(avg_lengthFluc,[(16:-2:2),15:-2:1],error_lengthFluc,'k');
% xticklabel_rotate(1:length(avg_fluc_per_nuc),-45,mutant_cat1,'interpreter','none')
set(gca,'YTick',1:16,'YTicklabel',nameIntercleave(end:-1:1));
title('average sumed length of fluctuations detected per nuclei (errorbar standard_error)')
xlabel('average total fluctutaion length per nucleus')
% print(gcf,[savepicdir,'\sumedfluclength_per_nuc'],'-deps');
print(gcf,[savepicdir,'\sumedfluclength_per_nuc'],'-dpng');

% average peak value of all fluctuations
figure(1003)
barh(18-2*(1:numM),avg_maxheights(1:8),0.5,'r'); hold on;
barh(17-2*(1:numM),avg_maxheights(9:16),0.5,'b');hold on;
herrorbar(avg_maxheights,[(16:-2:2),15:-2:1],error_maxheights,'k');
% xticklabel_rotate(1:length(avg_fluc_per_nuc),-45,mutant_cat1,'interpreter','none')
set(gca,'YTick',1:16,'YTicklabel',nameIntercleave(end:-1:1));
title('average of max height of fluctuations (errorbar std)')
xlabel('relative fluctuation size')
% print(gcf,[savepicdir,'\avgmaxheight'],'-deps');
print(gcf,[savepicdir,'\avgmaxheight'],'-dpng');

% average duration of all fluctuations
figure(1004)
barh(18-2*(1:numM),avg_durations(1:8),0.5,'r'); hold on;
barh(17-2*(1:numM),avg_durations(9:16),0.5,'b');hold on;
herrorbar(avg_durations,[(16:-2:2),15:-2:1],error_durations,'k');
% xticklabel_rotate(1:length(avg_fluc_per_nuc),-45,mutant_cat1,'interpreter','none')
set(gca,'YTick',1:16,'YTicklabel',nameIntercleave(end:-1:1));
title('average durations of fluctuations (errorbar std)')
xlabel('seconds')
% print(gcf,[savepicdir,'\avgdurations'],'-deps');
print(gcf,[savepicdir,'\avgdurations'],'-dpng');

% average mean width of all fluctuations
figure(1005)
barh(18-2*(1:numM),avg_meanwidth(1:8),0.5,'r'); hold on;
barh(17-2*(1:numM),avg_meanwidth(9:16),0.5,'b');hold on;
herrorbar(avg_meanwidth,[(16:-2:2),15:-2:1],error_meanwidth,'k');
% xticklabel_rotate(1:length(avg_fluc_per_nuc),-45,mutant_cat1,'interpreter','none')
set(gca,'YTick',1:16,'YTicklabel',nameIntercleave(end:-1:1));
title('average mean width of fluctuations (errorbar std)')
xlabel('degree')
% print(gcf,[savepicdir,'\avgmeanwidth'],'-deps');
print(gcf,[savepicdir,'\avgmeanwidth'],'-dpng');

tdata=table(nameAll',avg_numFluc',error_numFluc',...avg_lengthFluc',error_lengthFluc',...
    avg_durations',error_durations',avg_maxheights',error_maxheights',...
    avg_meanwidth',error_meanwidth',...
    'VariableNames',{'Name','average_number_of_fluctuations_per_nucleus',...
    'standard_error_of_number_of_fluctuations_per_nucleus',...
    'average_durations','standard_deviation_of_durations',...
    'average_maxinum_fluctuation_height','standard_deviation_of_max_fluctuation_height',...
    'average_mean_fluctuation_width','standard_devation_of_mean_fluctuation_width'});
write(tdata,[savepicdir,'\stats.xls']);

%% correlation scatter plot
% non MBC all maxheight vs durations
figure(2001);
set(gcf,'Position',[0 0 1400 800]);
for icat=1:numM
    scatter(durations{icat},maxheights{icat});hold on;
end
xlabel('durations of fluctuations(s)');
ylabel('max height of relative fluctuations');
axis([0 150 0 0.4])
legend(nameAll(1:8));
print(gcf,[savepicdir,'\maxheightvsduration_nonMBC'],'-dpng');
print(gcf,[savepicdir,'\maxheightvsduration_nonMBC'],'-deps');

% MBC all maxheight vs durations
figure(2002);
set(gcf,'Position',[0 0 1400 800]);
for icat=numM+1:2*numM
    scatter(durations{icat},maxheights{icat});hold on;
end
xlabel('durations of fluctuations(s)');
ylabel('max height of relative fluctuations');
axis([0 150 0 0.4])
legend(nameAll(9:16));
print(gcf,[savepicdir,'\maxheightvsduration_MBC'],'-dpng');
print(gcf,[savepicdir,'\maxheightvsduration_MBC'],'-deps');

% pair all maxheight vs durations
figure(2003);
set(gcf,'Position',[0 0 1400 800]);
for icat=1:numM
    clf;
    scatter(durations{icat},maxheights{icat},'rs');hold on;
    scatter(durations{icat+numM},maxheights{icat+numM},'bo');hold on;
    xlabel('durations of fluctuations(s)');
    ylabel('max height of relative fluctuations');
    axis([0 150 0 0.4])
    legend(nameAll([icat,icat+numM]));
    print(gcf,[savepicdir,'\maxheightvsduration_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\maxheightvsduration_',nameAll{icat}],'-deps');
end

%% correlation scatter plot
% non MBC all maxheight vs meanwidth
figure(2001);
set(gcf,'Position',[0 0 1400 800]);
for icat=1:numM
    scatter(meanwidth{icat},maxheights{icat});hold on;
end
xlabel('mean width of fluctuations(degree)');
ylabel('max height of relative fluctuations');
axis([10 90 0 0.4])
legend(nameAll(1:8));
print(gcf,[savepicdir,'\maxheightvsmw_nonMBC'],'-dpng');
print(gcf,[savepicdir,'\maxheightvsmw_nonMBC'],'-deps');

% MBC all maxheight vs meanwidth
figure(2002);
set(gcf,'Position',[0 0 1400 800]);
for icat=numM+1:2*numM
    scatter(meanwidth{icat},maxheights{icat});hold on;
end
xlabel('mean width of fluctuations(degree)');
ylabel('max height of relative fluctuations');
axis([10 90 0 0.4])
legend(nameAll(9:16));
print(gcf,[savepicdir,'\maxheightvsmw_MBC'],'-dpng');
print(gcf,[savepicdir,'\maxheightvsmw_MBC'],'-deps');

% pair all maxheight vs meanwidth
figure(2003);
set(gcf,'Position',[0 0 1400 800]);
for icat=1:numM
    clf;
    scatter(meanwidth{icat},maxheights{icat},'rs');hold on;
    scatter(meanwidth{icat+numM},maxheights{icat+numM},'bo');hold on;
    xlabel('mean width of fluctuations(degree)');
    ylabel('max height of relative fluctuations');
    axis([10 90 0 0.4])
    legend(nameAll([icat,icat+numM]));
    print(gcf,[savepicdir,'\maxheightvsmw_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\maxheightvsmw_',nameAll{icat}],'-deps');
end


%% calculate all heights and widths of fluctuation
% pkh_th=0.5;
max_tl=201;
mid_tl=(max_tl+1)/2;
% all heights
all_heights=cell(size(alltracks));
all_widths=cell(size(alltracks));
all_heights_norm=cell(size(alltracks));
% max_heights=cell(size(alltracks));
for icat=1:numM*2
    subtracks=alltracks{icat};
    height_save=cell(1,max_tl);
    width_save=cell(1,max_tl);
    height_norm_save=cell(1,max_tl);
    %     mh_save=zeros(1,length(subtracks));
    for ifluc=1:length(subtracks)
        heights=subtracks{ifluc}(:,4)';
        heights_norm=heights/max(heights);
        widths=subtracks{ifluc}(:,5)';
        ts=1:length(heights);
        [mh,maxheighttime]=max(heights);
        ts=ts-maxheighttime+mid_tl;
        %         mh_save(ifluc)=mh;
        for its=1:length(ts);
            height_save{ts(its)}=[height_save{ts(its)},heights(its)];
            height_norm_save{ts(its)}=[height_norm_save{ts(its)},heights_norm(its)];
            width_save{ts(its)}=[width_save{ts(its)},widths(its)];
        end
    end
    all_heights{icat}=height_save;
    all_widths{icat}=width_save;
    all_heights_norm{icat}=height_norm_save;
end
%average heights and standard deviation
avg_heights=zeros(max_tl,numM*2);
std_heights=zeros(max_tl,numM*2);
error_heights=zeros(max_tl,numM*2);
for icat=1:numM*2
    avg_heights(:,icat)=cellfun(@mean,all_heights{icat})';
    std_heights(:,icat)=cellfun(@std,all_heights{icat})';
    error_heights(:,icat)=std_heights(:,icat)./sqrt(cellfun(@length,all_heights{icat})');
end

%average normal heights and standard deviation
avg_heights_norm=zeros(max_tl,numM*2);
std_heights_norm=zeros(max_tl,numM*2);
error_heights_norm=zeros(max_tl,numM*2);
for icat=1:numM*2
    avg_heights_norm(:,icat)=cellfun(@mean,all_heights_norm{icat})';
    std_heights_norm(:,icat)=cellfun(@std,all_heights_norm{icat})';
    error_heights_norm(:,icat)=std_heights_norm(:,icat)./sqrt(cellfun(@length,all_heights_norm{icat})');
end

%average widths and standard deviation
avg_widths=zeros(max_tl,numM*2);
std_widths=zeros(max_tl,numM*2);
error_widths=zeros(max_tl,numM*2);
for icat=1:numM*2
    avg_widths(:,icat)=cellfun(@mean,all_widths{icat})';
    std_widths(:,icat)=cellfun(@std,all_widths{icat})';
    error_widths(:,icat)=std_widths(:,icat)./sqrt(cellfun(@length,all_widths{icat})');
end


%% fall rise time analysis, 

% fall rise of each mutant
figure(3001)
set(gcf,'Position',[0 0 1200 800])
for icat=1:numM*2
    clf
    avg_heights1=avg_heights(:,icat);
    error_heights1=error_heights(:,icat);
    ts=(1:mid_tl)*2.5;
    errorbar(ts,avg_heights1(mid_tl:end),error_heights1(mid_tl:end),'r');hold on;
    errorbar(ts,avg_heights1(mid_tl:-1:1),error_heights1(mid_tl:-1:1),'b');
    legend('fall','rise')
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
    title('average rise and fall curve');
    title(strrep(nameAll{icat},'_',' '));
    axis([1 100 -0.01 0.15])
    print(gcf,[savepicdir,'\risefall_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\risefall_',nameAll{icat}],'-deps');
end

% norm fall rise of each mutant
figure(3002)
set(gcf,'Position',[0 0 1200 800])
for icat=1:numM*2
    clf
    avg_heights1=avg_heights_norm(:,icat);
    error_heights1=error_heights_norm(:,icat);
    ts=(1:mid_tl)*2.5;
    errorbar(ts,avg_heights1(mid_tl:end),error_heights1(mid_tl:end),'r');hold on;
    errorbar(ts,avg_heights1(mid_tl:-1:1),error_heights1(mid_tl:-1:1),'b');
    legend('fall','rise')
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
    title('average rise and fall curve');
    title(strrep(nameAll{icat},'_',' '));
    axis([1 100 -0.1 1])
    print(gcf,[savepicdir,'\risefallnorm_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\risefallnorm_',nameAll{icat}],'-deps');
end

% norm fall rise of each mutant
figure(3003)
set(gcf,'Position',[0 0 1200 800])
for icat=1:numM
    clf
    avg_heights1=avg_heights_norm(:,icat);
    error_heights1=error_heights_norm(:,icat);
    avg_heights2=avg_heights_norm(:,icat+numM);
    error_heights2=error_heights_norm(:,icat+numM);
    ts=(1:mid_tl)*2.5;
    errorbar(ts,avg_heights1(mid_tl:end),error_heights1(mid_tl:end),'r');hold on;
    errorbar(ts,avg_heights1(mid_tl:-1:1),error_heights1(mid_tl:-1:1),'y');
    errorbar(ts,avg_heights2(mid_tl:end),error_heights2(mid_tl:end),'b');hold on;
    errorbar(ts,avg_heights2(mid_tl:-1:1),error_heights2(mid_tl:-1:1),'g');
    legend('fall','rise','MBC fall','MBC rise')
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
    title('average rise and fall curve');
    title(strrep(nameAll{icat},'_',' '));
    axis([1 100 -0.1 1])
    print(gcf,[savepicdir,'\2risefallnorm_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\2risefallnorm_',nameAll{icat}],'-deps');
end

% norm rise of each mutant
figure(3004)
set(gcf,'Position',[0 0 1200 800])
for icat=1:numM
    clf
    avg_heights1=avg_heights_norm(:,icat);
    error_heights1=error_heights_norm(:,icat);
    avg_heights2=avg_heights_norm(:,icat+numM);
    error_heights2=error_heights_norm(:,icat+numM);
    ts=(1:mid_tl)*2.5;
    errorbar(ts,avg_heights1(mid_tl:-1:1),error_heights1(mid_tl:-1:1),'y');hold on;
    errorbar(ts,avg_heights2(mid_tl:-1:1),error_heights2(mid_tl:-1:1),'g');hold on;
    legend('rise','MBC rise')
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
    title('average rise and fall curve');
    title(strrep(nameAll{icat},'_',' '));
    axis([1 100 -0.1 1])
    print(gcf,[savepicdir,'\2risenorm_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\2risenorm_',nameAll{icat}],'-deps');
end

% norm fall of each mutant
figure(3005)
set(gcf,'Position',[0 0 1200 800])
for icat=1:numM
    clf
    avg_heights1=avg_heights_norm(:,icat);
    error_heights1=error_heights_norm(:,icat);
    avg_heights2=avg_heights_norm(:,icat+numM);
    error_heights2=error_heights_norm(:,icat+numM);
    ts=(1:mid_tl)*2.5;
    errorbar(ts,avg_heights1(mid_tl:end),error_heights1(mid_tl:end),'r');hold on;
    errorbar(ts,avg_heights2(mid_tl:end),error_heights2(mid_tl:end),'b');hold on;
    legend('fall','MBC fall')
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
    title('average rise and fall curve');
    title(strrep(nameAll{icat},'_',' '));
    axis([1 100 -0.1 1])
    print(gcf,[savepicdir,'\2fallnorm_',nameAll{icat}],'-dpng');
    print(gcf,[savepicdir,'\2fallnorm_',nameAll{icat}],'-deps');
end
%% all fluctuation average aligned
colorSet = hsv(8);
ts=(-mid_tl+1:mid_tl-1)'*2.5;

figure(4003)
set(gcf,'Position',[0 0 1200 800])
for i=1:8
plot(ts,avg_heights(:,i),'color',colorSet(i,:),'LineWidth',2);hold on;
end
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
axis([-50 50 -0.02 0.2])
legend(strrep(nameAll(1:8),'_',' '));
title('average fluctuation aligned at maximun height');
print(gcf,[savepicdir,'\risefall_NMBC'],'-dpng');
print(gcf,[savepicdir,'\risefall_NMBC'],'-deps');

figure(4004)
set(gcf,'Position',[0 0 1200 800])
% ts=(-mid_tl+1:mid_tl-1)'*ones(1,numM)*2.5;
for i=9:16
plot(ts,avg_heights(:,i),'color',colorSet(i-8,:),'LineWidth',2);hold on;
end
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
axis([-50 50 -0.02 0.2])
legend(strrep(nameAll(9:16),'_',' '));
title('average fluctuation aligned at maximun height');
print(gcf,[savepicdir,'\risefall_MBC'],'-dpng');
print(gcf,[savepicdir,'\risefall_MBC'],'-deps');

figure(4005)
set(gcf,'Position',[0 0 1200 800])
% ts=(-mid_tl+1:mid_tl-1)'*ones(1,numM)*2.5;
for i=1:8
plot(ts,avg_heights_norm(:,i),'color',colorSet(i,:),'LineWidth',2);hold on;
end
    xlabel('time(seconds)');
    ylabel('relative fluctuation');
axis([-50 50 -0.2 1])
legend(strrep(nameAll(1:8),'_',' '));
title('average normalized fluctuation aligned at maximun height');
print(gcf,[savepicdir,'\risefallnorm_NMBC'],'-dpng');
print(gcf,[savepicdir,'\risefallnorm_NMBC'],'-deps');

figure(4006)
set(gcf,'Position',[0 0 1200 800])
% ts=(-mid_tl+1:mid_tl-1)'*ones(1,numM)*2.5;
for i=9:16
plot(ts,avg_heights_norm(:,i),'color',colorSet(i-8,:),'LineWidth',2);hold on;
end
xlabel('time(seconds)');
    ylabel('relative fluctuation');
axis([-50 50 -0.2 1])
legend(strrep(nameAll(9:16),'_',' '));
title('average normalized fluctuation aligned at maximun height');
print(gcf,[savepicdir,'\risefallnorm_MBC'],'-dpng');
print(gcf,[savepicdir,'\risefallnorm_MBC'],'-deps');


%% stats on rise
ts=(1:mid_tl)'*2.5;

fall_avg=array2table([ts,avg_heights_norm(mid_tl:end,:)]);
fall_avg.Properties.VariableNames=[{'timestep'},nameAll];
write(fall_avg,[savepicdir,'\fall_avg.xls']);

rise_avg=array2table([ts,avg_heights_norm(mid_tl:-1:1,:);]);
rise_avg.Properties.VariableNames=[{'timestep'},nameAll];
write(rise_avg,[savepicdir,'\rise_avg.xls']);

fall_error=array2table([ts,error_heights_norm(mid_tl:end,:)]);
fall_error.Properties.VariableNames=[{'timestep'},nameAll];
write(fall_error,[savepicdir,'\fall_error.xls']);

rise_error=array2table([ts,error_heights_norm(mid_tl:-1:1,:)]);
rise_error.Properties.VariableNames=[{'timestep'},nameAll];
write(rise_error,[savepicdir,'\rise_error.xls']);

%     'VariableNames',{'Name','average_number_of_fluctuations_per_nucleus',...
%     'standard_error_of_number_of_fluctuations_per_nucleus',...
%     'average_durations','standard_deviation_of_durations',...
%     'average_maxinum_fluctuation_height','standard_deviation_of_max_fluctuation_height',...
%     'average_mean_fluctuation_width','standard_devation_of_mean_fluctuation_width'});

%% pair comparison of all fluc trj

figure('Position',[0 0 1400 800])
for icat=1:numM
    clf
    for isub=1:2
        icat1=icat+(isub-1)*numM;
        subtracks=alltracks{icat1};
        subplot(1,2,isub);
        for ifluc=1:length(subtracks)
            heights=subtracks{ifluc}(:,4);
            ts=1:length(heights);
            [~,maxheighttime]=max(heights);
            tss=(ts-maxheighttime)*2.5;
            plot(tss,heights);    hold on;
        end
        %plot standard and mean
        allts=((1:max_tl)-mid_tl)'*ones(1,length(alltracks))*2.5;
        choosedata=icat1;
        allts1=allts(:,choosedata);
        avg_heights1=avg_heights(:,choosedata);
        error_heights1=error_heights(:,choosedata);
        errorbar(allts1,avg_heights1,error_heights1,'g');
        
        title(strrep(nameAll{icat1},'_',' '));
        axis([-80 80 -0.1 0.4]);
        xlabel('time(s)');
        ylabel('relative fluctuation size');
    end
    print(gcf,[savepicdir,'\heightpair_',nameAll{icat}],'-dpng');
end
% all trj in eps
figure('Position',[0 0 1400 800])
for icat=1:numM*2
    clf
        subtracks=alltracks{icat};
        for ifluc=1:length(subtracks)
            heights=subtracks{ifluc}(:,4);
            ts=1:length(heights);
            [~,maxheighttime]=max(heights);
            tss=(ts-maxheighttime)*2.5;
            plot(tss,heights);    hold on;
        end
        %plot standard and mean
        allts=((1:max_tl)-mid_tl)'*ones(1,length(alltracks))*2.5;
        errorbar(allts(:,icat),avg_heights(:,icat),error_heights(:,icat),'g');
        title(strrep(nameAll{icat},'_',' '));
        axis([-80 80 -0.1 0.4]);
        xlabel('time(s)');
        ylabel('relative fluctuation size');
    print(gcf,[savepicdir,'\alltrj_',nameAll{icat}],'-deps');
end
%% pair comparison of fluc widths 

figure('Position',[50 50 1200 800])
for icat=1:numM
    clf
    for isub=1:2
        icat1=icat+(isub-1)*numM;
        subtracks=alltracks{icat1};
        subplot(1,2,isub);
        for ifluc=1:length(subtracks)
            heights=subtracks{ifluc}(:,4);
            widths=subtracks{ifluc}(:,5);
            ts=1:length(heights);
            [~,maxheighttime]=max(heights);
            tss=(ts-maxheighttime)*2.5;
            validind=find(widths>0);
            plot(tss(validind),widths(validind));    hold on;
        end        
        title(strrep(nameAll{icat1},'_',' '));
        axis([-50 50 0 20]);
        xlabel('time(s)');
        ylabel('relative fluctuation size');
    end
    print(gcf,[savepicdir,'\widthpair_',nameAll{icat}],'-dpng');
end


%%
%     numFluc=cellfun(@length,alltracks);
%     sumlengthFluc=cellfun(@(x)sum(cellfun(@(y)sum(y(:,5)~=0),x)),alltracks);
%     numNuc=zeros(size(numFluc));
%
% avg_fluc_per_nuc=numFluc./numNuc;
% avg_sumlength_fluc=sumlengthFluc./numNuc;
% afpnNMBC=avg_fluc_per_nuc(1:8);
% afpnMBC=avg_fluc_per_nuc(9:16);
% attlfNMBC=avg_sumlength_fluc(1:8);
% attlfMBC=avg_sumlength_fluc(9:16);
%
% % mutant name
% figure(1)
% barh(18-2*(1:numM),afpnNMBC,0.5,'r'); hold on;
% barh(17-2*(1:numM),afpnMBC,0.5,'b');
% set(gca,'Ytick',1:length(avg_fluc_per_nuc));
% % xticklabel_rotate(1:length(avg_fluc_per_nuc),-45,mutant_cat1,'interpreter','none')
% set(gca,'YTick',1:16,'YTicklabel',nameIntercleave(end:-1:1));
% title('average number of fluctuations detected per nuclei')
% xlabel('average number')
% print(gcf,[path,'\pic\fluc_per_nuc'],'-dpng');
%
% figure(2)
% barh(18-2*(1:numM),attlfNMBC,0.5,'r'); hold on;
% barh(17-2*(1:numM),attlfMBC,0.5,'b');
% set(gca,'Ytick',1:length(avg_fluc_per_nuc));
% % xticklabel_rotate(1:length(avg_fluc_per_nuc),-45,mutant_cat1,'interpreter','none')
% set(gca,'YTicklabel',nameIntercleave(end:-1:1));
% title('average sumed length of fluctuations detected per nuclei')
% xlabel('average number')
% print(gcf,[path,'\pic\sumedfluclength_per_nuc'],'-dpng');


%% plot all the fluctuation
% f=figure('Position',[50 50 900 600]);
% nameNMBC=mutant_cat(cellfun(@isempty,strfind(mutant_cat,'MBC')));
% nameMBC=mutant_cat(cellfun(@(x)~isempty(x),strfind(mutant_cat,'MBC')));
% trjNMBC=alltracks0(cellfun(@isempty,strfind(mutant_cat,'MBC')));
% trjMBC=alltracks0(cellfun(@(x)~isempty(x),strfind(mutant_cat,'MBC')));
% numM=length(nameMBC);
% wspace=0.01;
% for icat=1:length(nameMBC)
%     subtracks=trjNMBC{icat};
%     axes('Position',[wspace (icat-1)/numM+wspace 1/2-2*wspace 1/numM-2*wspace])
%     for ifluc=1:length(subtracks)
%         heights=subtracks{ifluc}(:,4);
%         ts=1:length(heights);
%         [~,maxheighttime]=max(heights);
%         ts=ts-maxheighttime;
%         plot(ts,heights);    hold on;
%         axis([-20 20 -0.1 0.4]);
%     end
%
%     subtracks=trjMBC{icat};
%     axes('Position',[1/2+wspace (icat-1)/numM+wspace 1/2-2*wspace 1/numM-2*wspace])
%     for ifluc=1:length(subtracks)
%         heights=subtracks{ifluc}(:,4);
%         ts=1:length(heights);
%         [~,maxheighttime]=max(heights);
%         ts=ts-maxheighttime;
%         plot(ts,heights);    hold on;
%         axis([-20 20 -0.1 0.4]);
%     end
%
% end
