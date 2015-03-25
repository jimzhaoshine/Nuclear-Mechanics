%spherical fit, use spherical fit to determine the center of the nuclei
% result turned out to be quite noisy
points=nm.points;
faces=nm.faces;
z_th=0.7;
in_center=(points(:,3)<z_th & points(:,3)> -z_th);
norm_c=points((in_center),:);

faces_c=[];
for i=1:size(faces,1)
    if sum(in_center(faces(i,:)))==3
        faces_c=[faces_c;faces(i,:)];
    end
end
%
% p.vertices=points;
% p.faces=faces_c;
% patch(p);

%% get the x y z ratio
zxr=nm.vox/nm.pix*1.33/1.516;
for inuc=1%:nm.num_nuc
    origin=[];
    sph_center=[];
    for iframe=1:nm.endframe
        nuc=nm.nuclei{iframe,inuc};
        %% 3d spherical fit
%         z_aberation_ratio=1.33/1.56;
%         pts=norm_c.*(nuc.r(in_center)*[1 1 z_aberation_ratio]);
%         np=size(pts,1);
%         ini3=[0 0 0 10];
%         lb3=[-2 -2 -2 1];
%         ub3=[2 2 2 100];
%         options_sphere=optimset('TolX',5e-2,'TolFun',1e-2,'Display','off');
%         sphere_error=@(P) sqrt((pts(:,1)-P(1)).^2+(pts(:,2)-P(2)).^2+...
%             (pts(:,3)-P(3)).^2)-P(4);
%         sphfit=lsqnonlin(sphere_error,ini3,lb3,ub3,options_sphere);
%         origin=[origin;nuc.origin];
%         sph_center=[sph_center;nuc.origin+[sphfit([1,2]) sphfit(3)/z_aberation_ratio]];

%         clf
%         patch1.vertices=points.*(nuc.r*[1 1 z_aberation_ratio]);
%         patch1.faces=faces_c;
%         patch2.vertices=points*sphfit(4)+ones(size(points,1),1)*sphfit(1:3);
%         patch2.faces=faces_c;
%         patch(patch1,'FaceColor','red','EdgeColor','none','FaceAlpha',.3);hold on;
%         patch(patch2,'FaceColor','blue','EdgeColor','none','FaceAlpha',.3);hold on;
%         view([1 1 0]);
%         axis([-10 10 -10 10 -10 10])
%         daspect([1 1 1])
%         grid off
%         axis off
%         camlight
%         lighting gouraud
%         pause
        %% area fit
        area=zeros(1,nm.sizeZ);
        for i=1:length(area)
            area(i)=nuc.contour(i).area;
        end
        stack_fit=(-2:2)+round(nuc.origin(3));
        stack_fit=max(1,stack_fit(1)):min(10,stack_fit(end));
        z=stack_fit*zxr;
        area_fit=area(stack_fit);
        ini_g=[10,nuc.origin(3)*zxr,1];
        lb_g=[1,1*zxr,0.1];
        ub_g=[30,10*zxr,10];
        options_sphere=optimset('TolX',5e-2,'TolFun',1e-2,'Display','off');
        gaussfit=@(P) sphere_area(P,z) - area_fit;
        gfit=lsqnonlin(gaussfit,ini_g,lb_g,ub_g,options_sphere);
        z_g=(1:.1:nm.sizeZ)*zxr;
        area_g=sphere_area(gfit,z_g);
        clf
        plot((1:length(area))*zxr,area,'ob',z,area_fit,'og',...
            z_g,area_g,'r-');
        gfit(3)
        pause
    end
    %% plot center
%     plot3(origin(:,1),origin(:,2),origin(:,3),'r-.');hold on;
%     plot3(sph_center(:,1),sph_center(:,2),sph_center(:,3),'g-.');hold on;
%     plot(origin-ones(size(origin,1),1)*mean(origin));hold on;
%     plot(sph_center-ones(size(sph_center,1),1)*mean(sph_center),'--');hold on;
    % conclusion:
    % using spherical fit bring in more noise then just the centroid
end

%%
for inuc=1
    for iframe=1
        nuc=nm.nuclei{iframe,inuc};
        pts=points.*(nuc.r*[1 1 1]);
        ptsc=pts(in_center,:);
    end
end