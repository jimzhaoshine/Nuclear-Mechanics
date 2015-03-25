%% refresh matlab
close all;
clear all;
clc;

% set location of dv movies
rootpath='C:\nuclei\data\';
% rootpath='E:\nuclei';
strainpath='wild type MBC';
moviename='sp10_MBC_05.mat';
% nm=nucmem3(fullfile([rootpath,strainpath],moviename));
load(fullfile([rootpath,strainpath],moviename));
nm.loadmovie(101*10);
%%
nm.remove_badcentroid(1);
nm.remove_badcentroid(51);
nm.remove_badcentroid(101);
nm.initialize;
%%
nm.pickorientation;
%%
nm.remove_badcontour(1);
%%
nm.process_allframes;
%% display result
nm.display_contour
%% remove contour
nm.remove_badcontour(1);
%% test zone
% nm.process_singleframe(1);
%% open file folder
winopen(nm.path)




