%% refresh matlab
close all;
clear all;
clc;

% set location of dv movies
rootpath='C:\nuclei\data\';
% rootpath='E:\nuclei';
strainpath='wild type MBC';
moviename='sp10_MBC_08.dv';

nm=nucmem3(fullfile([rootpath,strainpath],moviename));
% load movie
nm.endframe=101;
nm.loadmovie(101*10);
% select the threshold for analysis
 nm.get_centroid_firstframe;

%  % loadcentroid(nm);
% % choose which frame to end the analysis, default=100;
% nm.choose_endframe;
% % remove off focus, drifting, mitotic nucleus
nm.remove_badcentroid(1);
nm.remove_badcentroid(50);
nm.remove_badcentroid(101);
% initialze
nm.initialize();
%% unload movie and save
nm.save_contour;
%%
nm.pickorientation;
%%
nm.remove_badcontour(1);
    nm.process_allframes;
    nm.correct_drift(0.2,1);
    %%    
    for iframe=1:101
        display(['processing frame ',num2str(iframe),' of ',num2str(nm.endframe)]);
        nm.process_singleframe2(iframe);
    end
    nm.centralband_all;
    nm.save_contour(1);

%% display result
nm.display_contour
%% remove contour
nm.remove_badcontour(1);
%% test zone
% nm.process_singleframe(1);
%% open file folder
winopen(nm.path)




