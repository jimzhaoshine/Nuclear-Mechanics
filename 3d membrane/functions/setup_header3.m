% rootpath='C:\nuclei\post analysis result_0.2';
namegroup1={'wildtype','heh1','heh2','ima1'};
namegroup2={'wildtype_TSA','heh1_TSA','heh2_TSA','ima1_TSA'};
nameIntercleave=reshape([namegroup1;namegroup2],1,2*length(namegroup1));
nameAll=[namegroup1,namegroup2];
nameAll2=cellfun(@(x)strrep(x,'_',' '),nameAll,'UniformOutput',0);
colorAll=[255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162;255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162]/255;
numM=length(namegroup1);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);