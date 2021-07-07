function [ h2,b2 ] = parasort( h1,b1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
b2=b1;h2=h1;
temp=0;
for hindex1=1:length(h2)
    for hindex2=2:length(h2)
        if (h2(hindex2)<=h2(hindex2-1))
            temp=h2(hindex2);
            h2(hindex2)=h2(hindex2-1);
            h2(hindex2-1)=temp;
            temp=b2(hindex2);
            b2(hindex2)=b2(hindex2-1);
            b2(hindex2-1)=temp;
        end
    end
end

end


