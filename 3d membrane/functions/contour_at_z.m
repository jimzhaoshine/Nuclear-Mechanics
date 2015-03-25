function [ contour ] = contour_at_z( interp_z,nuc,points,edges)
%interpolate at z position

    x=nuc.r.*points(:,1)+nuc.origin(1);
    y=nuc.r.*points(:,2)+nuc.origin(2);
    z=nuc.r.*points(:,3)+nuc.origin(3);
    
    p1=edges(:,1);
    p2=edges(:,2);
    
    contour=[];
        cross_edge=edges((z(p1)<interp_z & z(p2)>=interp_z)|(z(p1)>interp_z & z(p2)<=interp_z),:);
        z1=z(cross_edge(:,1))-interp_z;
        z2=interp_z-z(cross_edge(:,2));
        w1=z2./(z1+z2);
        w2=z1./(z1+z2);
        xi=w1.*x(cross_edge(:,1))+w2.*x(cross_edge(:,2));
        yi=w1.*y(cross_edge(:,1))+w2.*y(cross_edge(:,2));
        cx=mean(xi);
        cy=mean(yi);
        xi=xi-cx;
        yi=yi-cy;
        [theta,rho]=cart2pol(xi,yi);
        res=sortrows([theta,rho],1);
        [xi,yi]=pol2cart(res(:,1),res(:,2));
        xi=xi+cx;
        yi=yi+cy;
        if ~isempty(xi)
            xi=[xi;xi(1)];
            yi=[yi;yi(1)];
        end
        contour.x=xi;
        contour.y=yi;
        contourS.area=polyarea(xi,yi);
end

