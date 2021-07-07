clc;clear all;
%hold on;
%%
f=fopen('SourceFFT.txt');
File=1;
l=fgetl(f);
in=1;
Prev=1600*(File-1);
while ischar(l)
    text{in}=l;
    data{in}=sscanf(text{in},'%f %f %f');
    fr(in)=data{in}(1);
    FFT(in)=data{in}(2)+1i*data{in}(3);
    
    l=fgetl(f);
    in=in+1;
end
figure
hold on;
plot(fr*(3e10),real(FFT))
plot(fr*(3e10),imag(FFT),'r')
plot(fr*(3e10),abs(FFT),'k')
axis([0 6e10 -1 1]);
legend('Real G(jw)','Imag G(jw)','|G(jw)|')
title('Fourier Transform of Gaussian Pulse')
xlabel('Frequency (Hz)')
