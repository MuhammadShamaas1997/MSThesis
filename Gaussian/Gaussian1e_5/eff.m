function [ out ] = eff( ax,ay )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for ind=1:length(ax)
%     ang0=atan(real(ay(ind))/real(ax(ind)));
%     abs0=hypot(real(ax(ind)),real(ay(ind)));
    out(ind)=real(ax(ind))+1i*real(ay(ind));
end
end

