points=nm.points;
faces=nm.faces;
%%
for inuc=1%:nm.num_nuc
    sum_r=0;
    for iframe=1:nm.endframe
        sum_r=sum_r+nm.nuclei{iframe,inuc}.r;
    end
end
mean_r=sum_r/nm.endframe;
mean_pts=[mean_r.*points(:,1),mean_r.*points(:,2),mean_r.*points(:,3)];
patchm.vertices=mean_pts;
patchm.faces=faces;
%%

for inuc=1%:nm.num_nuc
    figure 
    for iframe=1:nm.endframe
        clf
            nuc=nm.nuclei{iframe,inuc};
            pts=[nuc.r.*points(:,1),nuc.r.*points(:,2),nuc.r.*points(:,3)];
            patch1.vertices=pts;
            patch1.faces=faces;
            
            patch(patch1,'FaceColor','red','EdgeColor','none');hold on;
            patch(patchm,'FaceColor','blue','EdgeColor','none');hold on;
            view(3);
            axis([-10 10 -10 10 -10 10])
            daspect([1 1 1])
            grid off
            camlight
            lighting gouraud
        pause;
    end
end