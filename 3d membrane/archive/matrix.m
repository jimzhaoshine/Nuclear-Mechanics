[points,faces, edges,neighbors]=TriSphere(3);
num_points=size(points,1);
%%
arc_distance=zeros(num_points);
for j=1:num_points
    for i=1:num_points
        pi=points(i,:);
        pj=points(j,:);
        if i==j
            arc_distance(i,j)=0;
        elseif i>j
%             arc_distance(i,j)=sqrt(sum((pi-pj).^2));
            arc_distance(i,j)=acos(pi*pj');
        elseif i<j
            arc_distance(i,j)=arc_distance(j,i);
        end
    end
end
arc_distance=real(arc_distance);
precision=1e6;
arc_distance=round(arc_distance*precision)/precision;
SI(arc_distance);
%%


unique_arc=unique(arc_distance);
arc_ind=cell(1,length(unique_arc));
arc_num=zeros(1,length(unique_arc));
tracker=0;
for i=1:length(unique_arc)
    arc=unique_arc(i);
    [p1,p2]=find(arc_distance==arc);
    ind=[p1,p2];
    tracker=tracker+length(p1);
    ind=ind(p1<=p2);
    arc_ind{i}=ind;
    arc_num(i)=size(ind,1);
end
%%
plot(unique_arc,arc_num)

