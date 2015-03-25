function [ heights ] = risefall( P, times )
%1d generate height based on time
tmax=P(1);
t1=P(2);
t2=P(3);
hmax=P(4);
hbg=P(5);


heights=zeros(size(times));
for i=1:length(times)
    if times(i)<tmax-t1
        heights(i)=hbg;
    elseif times(i)<tmax
        heights(i)=hmax*(1-(tmax-times(i))/t1)+hbg;
    elseif times(i)<tmax+t2
        heights(i)=hmax*(1-(times(i)-tmax)/t2)+hbg;
    else
        heights(i)=hbg;
    end
end

end

