function [ img3,imgshift ] = shift_image( img, orientation, hpsiz )
%shift the image to the direction of the cell

img1=img(:,hpsiz+1:end-hpsiz);
if orientation>=0
    imgshift=round(orientation/360*64);
else
    imgshift=64+round(orientation/360*64);
end
img2=[img1(:,imgshift+1:end),img1(:,1:imgshift)];
img3=[img2(:,end-hpsiz+1:end),img2,img2(:,1:hpsiz)];

end

