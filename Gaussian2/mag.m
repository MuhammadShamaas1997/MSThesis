function [ amp ] = mag( a,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(a)
    amp(i)=sqrt(a(i)*a(i)+b(i)*b(i));
end
end

