%temp test for drift correcton
load('C:\nuclei\post analysis result_0.1\data\wild_type_04.mat');
inuc=3;
%% i should try 0.2 filter
fNorm=.2;
order = 5;      % order of filter -> lowest is 1, highest is whatever
[b1,a1] = butter(order, fNorm, 'low');
[b2,a2] = butter(order, 0.2, 'low');
zxr=nm.vox/nm.pix;
points=nm.points;
faces=nm.faces;

cnt=zeros(nm.endframe,3);
ocnt=zeros(nm.endframe,3);
for iframe=1:nm.endframe
    nuc=nm.nuclei{iframe,inuc};
    cnt(iframe,:)=nuc.center;
    ocnt(iframe,:)=nuc.origin;
end
fcnt =[ filtfilt(b1, a1, cnt(:,[1,2])),filtfilt(b2, a2, cnt(:,3))];
mcnt=mean(cnt,1);
mcnt=ones(size(cnt,1),1)*mcnt;
corcnt=fcnt-ocnt;
corcnt(:,3)=corcnt(:,3)*zxr;
plot((1:101)'*ones(1,6)*2.5,([cnt,fcnt]-[mcnt,mcnt])*p2um,'-');
legend('raw x','raw y','raw z','filtered x','filtered y','filtered z');
xlabel('time (s)');
ylabel('position (\mum)');
FigureFormat(gcf);
%%
cnt0=nm.nuclei{1,inuc}.origin;
f=figure('Position',[50 50 1400 700]);
framestep=5;
for iframe=1:nm.endframe-framestep
    %%
    nuc1=nm.nuclei{iframe,inuc};
    nuc2=nm.nuclei{iframe+framestep,inuc};
    clf
    pts1=[nuc1.r_new.*points(:,1),nuc1.r_new.*points(:,2),nuc1.r_new.*points(:,3)]...
        +ones(size(points,1),1)*nuc1.origin_new;
    pts1c=[nuc1.r_new.*points(:,1),nuc1.r_new.*points(:,2),nuc1.r_new.*points(:,3)];
    patch1.vertices=pts1*p2um;
    patch1.faces=faces;
    patch1c.vertices=pts1c*p2um;
    patch1c.faces=faces;
    pts2=[nuc2.r_new.*points(:,1),nuc2.r_new.*points(:,2),nuc2.r_new.*points(:,3)]...
        +ones(size(points,1),1)*nuc2.origin_new;
    pts2c=[nuc2.r_new.*points(:,1),nuc2.r_new.*points(:,2),nuc2.r_new.*points(:,3)];
    patch2.vertices=pts2*p2um;
    patch2.faces=faces;
    patch2c.vertices=pts2c*p2um;
    patch2c.faces=faces;
    
     axes('Position',[0.05 0.05 .45 .95]);
    patch(patch1,'FaceColor','g','EdgeColor','none','FaceAlpha',0.3);hold on;
    patch(patch2,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold on;
    view([0 0 1]);
    axis([cnt0(1)-10 cnt0(1)+10 cnt0(2)-10 cnt0(2)+10 cnt0(3)-10 cnt0(3)+10]*p2um)
    daspect([1 1 1])
    grid off
    camlight
    lighting gouraud
        legend('frame i','frame i+5')
        xlabel('x (um)')
        ylabel('y (um)');

     axes('Position',[.55 0.05 .45 .95]);
    patch(patch1c,'FaceColor','g','EdgeColor','none','FaceAlpha',0.3);hold on;
    patch(patch2c,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold on;
    view([0 0 1]);
    axis([-10 +10 -10 +10 -10 +10]*p2um)
    daspect([1 1 1])
    grid off
    camlight
    lighting gouraud
        legend('frame i','frame i+5')
        xlabel('x (um)')
        ylabel('y (um)');
 FigureFormat(gcf)
%     pause
end


