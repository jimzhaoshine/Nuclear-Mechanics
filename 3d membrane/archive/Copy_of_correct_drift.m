function [ nm ] = correct_drift( nm, varargin )

if nargin==1
    fNorm = .1;  % f/samplingf -> .01 very low pass filter  -> 1 no filter
elseif nargin==2
    fNorm=varargin{1};
end
order = 5;      % order of filter -> lowest is 1, highest is whatever
[b,a] = butter(order, fNorm, 'low');
zxr=nm.vox/nm.pix;
points=nm.points;
faces=nm.faces;

%% second neighbor list
neighbors=nm.neighbors;
neighbors2=nan(size(points,1),19);
for i=1:size(neighbors,1)
    nbtmp=[];
    for j=1:size(neighbors,2)
        if ~isnan(neighbors(i,j))
            nbtmp=[nbtmp,neighbors(neighbors(i,j),:)];
        end
    end
    nbtmp=unique(nbtmp);
    nbtmp=nbtmp(~isnan(nbtmp));
    neighbors2(i,1:length(nbtmp))=nbtmp;
end

px=[points(:,1);0];
py=[points(:,2);0];
pz=[points(:,3);0];
neighbors2(isnan(neighbors2))=length(px);
px2=px(neighbors2);
py2=py(neighbors2);
pz2=pz(neighbors2);
cosd=px(275)*px(277)+py(275)*py(277)+pz(275)*pz(277);
dtheta=sqrt(2*(1-cosd));

%% correct drift
for inuc=1:nm.num_nuc
    %% filter origin
    cnt=zeros(nm.endframe,3);
    ocnt=zeros(nm.endframe,3);
    for iframe=1:nm.endframe
        nuc=nm.nuclei{iframe,inuc};
        cnt(iframe,:)=nuc.center;
        ocnt(iframe,:)=nuc.origin;
    end
    fcnt = filtfilt(b, a, cnt);
    mcnt=mean(cnt,1);
    mcnt=ones(size(cnt,1),1)*mcnt;
    corcnt=fcnt-ocnt;
    corcnt(:,3)=corcnt(:,3)*zxr;
    %% interpolate
    for iframe=1:nm.endframe
        nuc=nm.nuclei{iframe,inuc};
        x=points(:,1).*nuc.r-corcnt(iframe,1);
        y=points(:,2).*nuc.r-corcnt(iframe,2);
        z=points(:,3).*nuc.r-corcnt(iframe,3);
        r=sqrt(x.^2+y.^2+z.^2);
        nx=(x./r)*ones(1,19);
        ny=(y./r)*ones(1,19);
        nz=(z./r)*ones(1,19);
        dist2=nx.*px2+ny.*py2+nz.*pz2;
        [sd2,sortI]=sort(dist2,2);
        i1=neighbors2(sub2ind(size(neighbors2),(1:length(x))',sortI(:,end-2)));
        i2=neighbors2(sub2ind(size(neighbors2),(1:length(x))',sortI(:,end-1)));
        i3=neighbors2(sub2ind(size(neighbors2),(1:length(x))',sortI(:,end-0)));
        w1=sqrt(2*(1-sd2(:,end-2)));
        w2=sqrt(2*(1-sd2(:,end-1)));
        w3=sqrt(2*(1-sd2(:,end-0)));
        ri=(w2.*w3.*r(i1)+w1.*w3.*r(i2)+w1.*w2.*r(i3))...
            ./(w2.*w3+w1.*w3+w1.*w2);
        fcenter=fcnt(iframe,:);
        nuc.r_new=ri;
        nuc.origin_new=fcenter;
        nm.nuclei{iframe,inuc}=nuc;
        %% test plot
        %         pts1=[nuc.r.*points(:,1),nuc.r.*points(:,2),nuc.r.*points(:,3)];
        %         patch1.vertices=pts1;
        %         patch1.faces=faces;
        %         pts2=[nuc.r_new.*points(:,1),nuc.r_new.*points(:,2),nuc.r_new.*points(:,3)];
        %         pts2=pts2+ones(size(pts2,1),1)*corcnt(iframe,:);
        %         patch2.vertices=pts2;
        %         patch2.faces=faces;
        %
        %         patch(patch1,'FaceColor','red','EdgeColor','none','FaceAlpha',.5);hold on;
        %         patch(patch2,'FaceColor','blue','EdgeColor','none','FaceAlpha',.5);hold on;
        %         view(3);
        %         axis([-10 10 -10 10 -10 10])
        %         daspect([1 1 1])
        %         grid off
        %         camlight
        %         lighting gouraud
    end
    clf
    plot(cnt-mcnt);
    hold on;
    plot(fcnt-mcnt);
    pause(1);
end

end