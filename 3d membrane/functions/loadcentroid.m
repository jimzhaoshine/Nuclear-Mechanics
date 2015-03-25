function [  ] = loadcentroid( nm )
%load centroid points from text file
[f1,f2,f3]=fileparts(nm.path);
filein=[f1,'\GoodNucleiIndex\',nm.filename,'.txt'];
M = dlmread(filein);
cnt_tmp1=[M(:,1),nm.sizeY-M(:,2),5.5+zeros(size(M,1),1)];

SI(nm.proj(1));hold on; plot(cnt_tmp1(:,1),cnt_tmp1(:,2),'.');

for i=1:size(cnt_tmp1,1)
    for j=1:size(nm.cnt_tmp,1)
        if sum((cnt_tmp1(i,:)-nm.cnt_tmp(j,:)).^2)<100
            cnt_tmp1(i,:)=nm.cnt_tmp(j,:);
        end
    end
end
nm.cnt_tmp=cnt_tmp1;

SI(nm.proj(1));hold on; plot(cnt_tmp1(:,1),cnt_tmp1(:,2),'.');

display(['centroid from file (',filein,') loaded.']);
end

