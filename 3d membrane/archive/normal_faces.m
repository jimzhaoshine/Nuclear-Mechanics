function [ n ] = normal_faces( points,faces )
%calculate normal values of faces
% points=nm.points;
% faces=nm.faces;
p1=points(faces(:,1),:);
p2=points(faces(:,2),:);
p3=points(faces(:,3),:);
edge1=p1-p2;
edge2=p2-p3;
n=cross(edge1,edge2);
n=n./(sqrt(sum(n.^2,2))*ones(1,3));

end

