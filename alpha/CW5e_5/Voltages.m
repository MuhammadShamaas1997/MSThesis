clc;clear all;
%hold on;
%%
f=fopen('Voltages.txt');
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
    
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f');
    dis=round(data{in}(2)*10+61);

    Im(data{in}(1),dis)=data{in}(3)+1i*data{in}(4);
    Vm(data{in}(1),dis)=data{in}(5)+1i*data{in}(6);
    
    l=fgetl(f);
    in=in+1;
end
