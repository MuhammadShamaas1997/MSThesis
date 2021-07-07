function [ out ] = runningmean( inm,num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
extra=((num-1)/2);
out=zeros(1,length(inm));
inm2=[zeros(1,extra) inm zeros(1,extra)];
for pl=(1+extra):(length(inm2)-extra)
sumn=mean(inm2((pl-extra):(pl+extra)));
    out(pl-extra)=sumn;
end

end

