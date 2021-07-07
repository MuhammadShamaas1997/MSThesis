clc;clear all;
%hold on;
%%
f=fopen('Energy.txt');
%f=fopen('Space.txt');
File=1;
l=fgetl(f);
in=1;
Prev=1600*(File-1);
while ischar(l)
%while (in <= (518*300))
%for kj=1:81
    %%disp(l);
    text{in}=l;
    
    data{in}=sscanf(text{in},'%f , %f , %f');
    
    Enx(in)=data{in}(1);
    Eny(in)=data{in}(2);
    Enz(in)=data{in}(3);
    
    l=fgetl(f);
    in=in+1;
end
