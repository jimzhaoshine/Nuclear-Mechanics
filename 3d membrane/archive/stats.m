%statistics of analyzed movie

rootpath='C:\nuclei';
folders=dir(rootpath);
categoryID=0;
for i=1:length(folders)
    
    strpath=folders(i).name;
%     if sum(strcmp(strpath,{'wild type' 'heh1' 'heh2' 'ima1'}))
        moviename=dir([rootpath,'\',strpath,'\*.mat']);
        if ~isempty(moviename)
            categoryID=categoryID+1;
            category=moviename(1).name;
            category=category(1:end-7);
            num_total=0;
            num_completed=0;
            nuc=[];
            for j=1:length(moviename)
                load([rootpath,'\',strpath,'\',moviename(j).name]);
                num_total=num_total+nm.num_nuc;
                if nm.continuefrom_frame>=102%nm.endframe
                    num_completed=num_completed+nm.num_nuc;
                else
                    moviename(j).name
                end
                nuc=[nuc,nm.nuclei];
            end
            report(categoryID).nuc=nuc;
            report(categoryID).name=category;
            report(categoryID).num_total=num_total;
            report(categoryID).num_completed=num_completed;
        end
%     end
end

%%
for catid=1:length(report)
nuc=report(catid).nuc;

area=zeros(size(nuc));
for j=1:size(nuc,2)
for i=1:size(nuc,1)
    area(i,j)=max([nuc{i,j}.contour.area]);
end
end
report(catid).area=area;

end


%%
barc=50:30:350;
area_count=[];
legend_str=[];
for catid=1:length(report)
    mean_area=mean(mean(report(catid).area));
    area_count(catid,:)=hist(mean(report(catid).area),barc);
    total_count=sum(area_count(catid,:));
    area_count(catid,:)=area_count(catid,:)/sum(area_count(catid,:));
    legend_str=[legend_str,{[report(catid).name,': ',num2str(mean_area),' ',num2str(total_count)]}];
end
bar(barc',area_count','hist')
legend(legend_str)
%%
plot(area)
