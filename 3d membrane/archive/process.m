allfile=dir('*.mat');
for i=1:length(allfile)
load(allfile(i).name);
nm.path(1)='C';
nm.loadmovie;

nm.process_singleframe(1);
for inuc=1:nm.num_nuc
    nm.cnt_tmp(inuc,:)=nm.nuclei{1,inuc}.center;
end
nm.process_singleframe(1);
for inuc=1:nm.num_nuc
    nm.cnt_tmp(inuc,:)=nm.nuclei{1,inuc}.center;
end
% process all frames
nm.process_allframes;

savefile=fullfile(nm.path,[nm.filename,'.mat']);
nm.mov=[];
save(savefile,'nm');

end