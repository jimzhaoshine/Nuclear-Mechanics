function [ nm ] = process_singleframe( nm ,iframe)
%process a single frame of the movie
showplot=1;
%% simulate a single frame
%initialize
num_nuc=nm.num_nuc;
currentcnt=nm.cnt_tmp;
wsz=nm.wsize;
hw=(wsz-1)/2;
rs=nm.rs;
points=nm.points;
faces=nm.faces;
edges=nm.edges;
neighbors=nm.neighbors;
zxr=nm.vox/nm.pix;
zexpand=nm.zexpand;

%get img
img=nm.grab3(iframe);
% img=bpass(img,.5,5);

for inuc=1:num_nuc
    tic
    %get center
    if iframe==1
    nuc_center=round(currentcnt(inuc,:));
    nuc_center(3)=round(currentcnt(inuc,3)-5.5)+5.5;
    else
        tmpcnt=nm.nuclei{iframe-1,inuc}.center;
        nuc_center(1:2)=round(tmpcnt(1:2));
        nuc_center(3)=round(tmpcnt(3)-5.5)+5.5;
    end
    %window image
    wimg=img(nuc_center(2)-hw:nuc_center(2)+hw,nuc_center(1)-hw:nuc_center(1)+hw,:);
    %interpolated intensity
    wimg1=[];
    zrange=(-4.5-zexpand:4.5+zexpand)+nuc_center(3);
    for iz=zrange
        if iz<1
            wimg1=cat(3,wimg1,wimg(:,:,1));
        elseif iz>10
            wimg1=cat(3,wimg1,wimg(:,:,10));
        else
            wimg1=cat(3,wimg1,wimg(:,:,iz));
        end
    end
    I=sum(wimg1(nm.linearindex).*nm.weights,3);
    
    %minimize energy
    rI=I';
    cost=100;
    rs0=rs(1)/(rs(2)-rs(1))-1;
    energy=@(ind)contour_energy_3(ind,rI,cost,rs0,edges);
    if iframe==1
    [~,ini_c]=max(rI);
    else
        ini_c=(nm.nuclei{iframe-1,inuc}.r)';
    end
    lb_c=zeros(1,size(points,1))+1;
    ub_c=zeros(1,size(points,1))+length(rs)-0.001;
    options=optimoptions('fmincon','MaxFunEvals', 1e6,'TolX',.01,'Display','final');
    [indr,fval,exitflag,output] = fmincon(energy,ini_c,[],[],[],[],lb_c,ub_c,[],options);
    [~,intensity]=energy(indr);
    r=rs(1)+(indr-1)*(rs(2)-rs(1));
    %update center
    x=r'.*points(:,1);
    y=r'.*points(:,2);
    z=r'.*points(:,3);
    meandcnt=[mean(x),mean(y),mean(z)/zxr];
%     currentcnt(inuc,:)=nuc_center+meandcnt(1:3);
    
    %calculate volume
    volume=trisphere_volume([x,y,z],faces);
    
    % calculate contour in each plane
    p1=edges(:,1);
    p2=edges(:,2);
    zstack=((1:10)-nuc_center(3))*zxr;
    contour=[];
    for i=1:nm.numstacks
        istack=zstack(i);
        cross_edge=edges((z(p1)<istack & z(p2)>=istack)|(z(p1)>istack & z(p2)<=istack),:);
        z1=z(cross_edge(:,1))-istack; 
        z2=istack-z(cross_edge(:,2));
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
        contour(i).x=xi;
        contour(i).y=yi;
        contour(i).area=polyarea(xi,yi);
    end

    
    % save data
    nuc.r=r';
    nuc.window_center=nuc_center;
    nuc.center=meandcnt+nuc_center;
    nuc.volume=volume;
    nuc.intensity=intensity;
    nuc.contour=contour;
    nm.nuclei{iframe,inuc}=nuc;
    
    toc
    %% plot 3d
    if showplot
        figure(101)
        nuc=nm.nuclei{iframe,inuc};
        pts=[nuc.r.*points(:,1),nuc.r.*points(:,2),nuc.r.*points(:,3)];
        TR = triangulation(faces,pts);
        trisurf(TR,'FaceColor','red','EdgeColor','black');
        view(3);
        daspect([1 1 1])
        camlight
        lighting gouraud
    end
    
    %% plot energy map compare
    if showplot
        figure(102)
        nuc=nm.nuclei{iframe,inuc};
        intensity=nuc.intensity;
        neighbors=nm.neighbors;
        neighbors(1:12,6)=(1:12)';
        r_energy=sum((indr(neighbors)-indr'*ones(1,6)).^2,2)*cost;
        plot(1:length(intensity),intensity,1:length(intensity),r_energy)
        legend('intensity','bending energy')
    end
    %% plot stack image
    if showplot
        f=figure(103);
        nuc=nm.nuclei{iframe,inuc};
        set(f,'Position',[50 50 1500 600]);
        for i=1:10
            xi=nuc.contour(i).x;
            yi=nuc.contour(i).y;
            pr=floor((i-1)/5);
            pc=i-5*pr-1;
            pr=1-pr;
            axes('Unit','pixel','Position',[pc*300 pr*300 300 300]);
            SI(wimg(:,:,i));
            hold on;
            plot(xi+hw+1,yi+hw+1,'-');
            axis off;
            box on;
            text(5,3,['Zstack:',num2str(i)],'Color','r','FontWeight','bold','FontSize',15);
            text(5,wsz-3,['Area:',num2str(contour(i).area)],'Color','y','FontWeight','bold','FontSize',15);
            if i==1
            text(5,10,[{'window center'},{['x:',num2str(nuc_center(1))]},...
                {['y:',num2str(nuc_center(2))]},{['z:',num2str(nuc_center(3))]}],...
                'Color','m','FontWeight','bold','FontSize',15);
            text(5,20,[{'nuc center'},{['x:',num2str(nuc.center(1))]},...
                {['y:',num2str(nuc.center(2))]},{['z:',num2str(nuc.center(3))]}],...
                'Color','g','FontWeight','bold','FontSize',15);
            end
        end
    end
end
%%

    
end


