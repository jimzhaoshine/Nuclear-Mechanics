

function singleworker(name)
    nmdata=load(name);
    nm=nmdata.nm;
    nm.loadmovie(101*10);
    nm.initialize;
    nm.process_allframes;
    nm.correct_drift(0.2,1);
    for iframe=1:101
        display(['processing(2nd) frame ',num2str(iframe),' of ',num2str(nm.endframe)]);
        nm.process_singleframe2(iframe);
    end
    nm.centralband_all;
    nm.detect_fluctuation;
    nm.save_contour(1);
end

