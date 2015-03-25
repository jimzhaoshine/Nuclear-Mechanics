% rootpath='C:\nuclei\post analysis result_0.2';
nameNMBC={'wild_type','clr4','swi6'};
nameMBC={'sp10_MBC','clr4_MBC','swi6_MBC'};
nameIntercleave=reshape([nameNMBC;nameMBC],1,2*length(nameNMBC));
nameAll=[nameNMBC,nameMBC];
nameAll2=cellfun(@(x)strrep(x,'_',' '),nameAll,'UniformOutput',0);
colorAll=[255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162;255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162]/255;
numM=length(nameNMBC);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);