function nuc=centralband(nuc,nm,points2,imsize,aberation)
%% convert 3D data to a 2D central band
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! no need to correct aberation if
%corrected in movie!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% interpolate to 2d image
% face centers

points=nm.points;
faces=nm.faces;

fc=(points(faces(:,1),:)+points(faces(:,2),:)+points(faces(:,3),:))/3;
fx2=ones(size(points2,1),1)*fc(:,1)';
fy2=ones(size(points2,1),1)*fc(:,2)';
fz2=ones(size(points2,1),1)*fc(:,3)';

%interpolate
r=nuc.r_new;
% norm value of each pixel
nx=points2(:,1)*ones(1,size(faces,1));
ny=points2(:,2)*ones(1,size(faces,1));
nz=points2(:,3)*ones(1,size(faces,1));
% find pixel in which triangle
dist2=sqrt((nx-fx2).^2+(ny-fy2).^2+(nz-fz2).^2);
[sd2,sortI]=sort(dist2,2);
i1=faces(sortI(:,1),1);
i2=faces(sortI(:,1),2);
i3=faces(sortI(:,1),3);
% old points
points_old=points.*(r*[1 1 aberation]);
p1=points_old(i1,:);
p2=points_old(i2,:);
p3=points_old(i3,:);
% intercept
facevec=cross(p3-p2,p2-p1,2);
facenorm=facevec./sqrt(sum(facevec.^2,2)*[1 1 1]);
d1=(sum(p1.*facenorm,2));
d2=(sum(points2.*facenorm,2));
pointsi=((d1./d2)*[1 1 1]).*points2;
ri=sqrt(sum(pointsi.^2,2));
nuc.img=reshape(ri,imsize);

% plot
if 0
    clf
    pts1=[ri.*points2(:,1),ri.*points2(:,2),ri.*points2(:,3)];
    patch1.vertices=pts1;
    patch1.faces=faces2;
    pts2=[nuc.r_new.*points(:,1),nuc.r_new.*points(:,2),nuc.r_new.*points(:,3)];
    %     pts2=pts2+ones(size(pts2,1),1)*corcnt(iframe,:);
    patch2.vertices=pts2;
    patch2.faces=faces;
    
    patch(patch1,'FaceColor','r','EdgeColor','none','FaceAlpha',.3);hold on;
    patch(patch2,'FaceColor','g','EdgeColor','none','FaceAlpha',.3);hold on;
    view(3);
    axis([-10 10 -10 10 -10 10])
    daspect([1 1 1])
    grid off
    camlight
    lighting gouraud
    pause
end

end


