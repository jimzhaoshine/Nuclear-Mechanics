% fourier analysis on each frames

% result : fourier analysis is not a good one identifying peak , should try
% paticle tracking more, and rembmer to take into the periodic condition of
% theta coordinate in fft

points=nm.points;
faces=nm.faces;
for inuc=1
    %% calculate mean contour
    sum_r=0;
    for iframe=1:nm.endframe
        sum_r=sum_r+nm.nuclei{iframe,inuc}.r_new;
    end
    mean_r=sum_r/nm.endframe;
    mean_pts=[mean_r.*points(:,1),mean_r.*points(:,2),mean_r.*points(:,3)];
    patchm.vertices=mean_pts;
    patchm.faces=faces;
    meannuc=nm.nuclei{1,inuc};
    meannuc.r_new=mean_r;
    [ points2,faces2,psi,theta ] = setup_img(  );
    meannuc=centralband(meannuc,nm,points2,size(psi));
    mean_img=meannuc.img;
    spt_fft_mean=zeros(19,64);
    for iframe=1:nm.endframe
        nuc=nm.nuclei{iframe,inuc};
        spt_fft=fft2(nuc.img-mean_img);
        spt_fft_mean=spt_fft_mean+abs(spt_fft);
    end
    spt_fft_mean=spt_fft_mean/nm.endframe;
%     spt_fft_mean(1)=0;
    SI(spt_fft_mean);
end
