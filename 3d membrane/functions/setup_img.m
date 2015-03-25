function [ points2,faces2,psi2,theta2 ] = setup_img(  )
%setup image area for 3d to 2d translation
num_points=64;
angle_step=2*pi/num_points;
psi1=(0:num_points-1)*angle_step;
theta1=(-num_points/7:num_points/7)*angle_step;
[psi2,theta2]=meshgrid(psi1,theta1);

%setup 3d angle for each pixel
x=reshape(cos(psi2).*cos(theta2),numel(psi2),1);
y=reshape(sin(psi2).*cos(theta2),numel(psi2),1);
z=reshape(sin(theta2),numel(psi2),1);
points2=[x,y,z];
faces2=[];
for i=2:size(psi2,1)
    for j=2:size(psi2,2)
        faces2=[faces2;sub2ind(size(psi2),[i-1 i-1 i i],[j-1 j j j-1])];
    end
end
end

