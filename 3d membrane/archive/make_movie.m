% function [  ] = make_movie( nm, inuc )

% daObj=VideoWriter('m1.avi','MPEG-4');
% daObj.set('Quality',100);
% daObj.FrameRate=1;
% open(daObj);

%make movie
inuc=1;
hw=(nm.wsize-1)/2;
points=nm.points;
f=figure;
for iframe=1:nm.endframe
    clf
    nuc=nm.nuclei{iframe,inuc};
    img=nm.grab3(iframe);
    round_origin=round(nuc.origin(1:2));
    d_origin=nuc.origin(1:2)-round_origin;
    wimg=img(round_origin(2)-hw:round_origin(2)+hw,...
        round_origin(1)-hw:round_origin(1)+hw,:);
    set(f,'Position',[0 50 1000 750]);
    
    axes('Unit','pixel','Position',[0 400 350 350]);
    pts=[nuc.r.*points(:,1),nuc.r.*points(:,2),nuc.r.*points(:,3)];
    TR = triangulation(nm.faces,pts);
    trisurf(TR,'FaceColor','red','EdgeColor','black');
    axis([-15 15 -15 15 -15 15]);
    grid off
    view([1 0 0]);
    daspect([1 1 1])
    camlight
    lighting gouraud
    
    axes('Unit','pixel','Position',[370 420 310 310]);
    intensity=nuc.intensity;
    neighbors=nm.neighbors;
    neighbors(1:12,6)=(1:12)';
    r_energy=sum((indr(neighbors)-indr'*ones(1,6)).^2,2)*cost;
    plot(1:length(intensity),intensity,1:length(intensity),r_energy)
    axis([0 length(intensity) 0 65535])
    legend('intensity','bending energy')
    
    axes('Unit','pixel','Position',[720 420 270 310]);
    area=zeros(1,nm.sizeZ);
    for i=1:length(area)
        area(i)=nuc.contour(i).area;
    end
    stack_fit=(-2:2)+round(nuc.origin(3));
    stack_fit=max(1,stack_fit(1)):min(10,stack_fit(end)); 
    area_fit=area(stack_fit);
    ini_g=[3,nuc.origin(3)];
    lb_g=[1,min(stack_fit)];
    ub_g=[10,max(stack_fit)];
    options_sphere=optimset('TolX',5e-2,'TolFun',1e-5,'Display','final');
    gaussfit=@(P) sphere_area(P,stack_fit) - area_fit;
    gfit=lsqnonlin(gaussfit,ini_g,lb_g,ub_g,options_sphere);
    stack_g=1:.1:nm.sizeZ;
    area_g=sphere_area(gfit,stack_g);
    plot(1:length(area),area,'ob',stack_g,area_g,'r-');
    axis([0 10 0 400]);
%     xlabel('stack id');
%     ylabel('area (pixel^2)');
    
    for i=1:10
        xi=nuc.contour(i).x;
        yi=nuc.contour(i).y;
        pr=floor((i-1)/5);
        pc=i-5*pr-1;
        pr=1-pr;
        imgsize=200;
        axes('Unit','pixel','Position',[pc*imgsize pr*imgsize imgsize imgsize]);
        SI(wimg(:,:,i));
        hold on;
        plot(xi+hw+1+d_origin(1),yi+hw+1+d_origin(2),'-','Linewidth',2);
        axis off;
        box on;
        text(5,3,['Zstack:',num2str(i)],'Color','r','FontWeight','bold','FontSize',15);
        text(5,wsz-3,['Area:',num2str(contour(i).area)],'Color','y','FontWeight','bold','FontSize',15);
        if i==1
            text(5,20,[{'nuc center'},{['x:',num2str(nuc.center(1))]},...
                {['y:',num2str(nuc.center(2))]},{['z:',num2str(nuc.center(3))]}],...
                'Color','g','FontWeight','bold','FontSize',15);
        end
    end
    pause(.2)
%     writeVideo(daObj,getframe(gcf));
end
close all;
% close(daObj);


% end

