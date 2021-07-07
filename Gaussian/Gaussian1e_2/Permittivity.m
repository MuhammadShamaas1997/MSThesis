clc;clear all;

fi=fopen('Permeability.txt');
l=fgetl(fi);
in=1;
while ischar(l)
%for kj=1:81
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f %f');
    f(in)=(data{in}(1));
    eps(in)=data{in}(2);
    
    l=fgetl(fi);
    in=in+1;
end

epsr=0.9999;
a0=1e-2;%0.1mm
c0=2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%3000GHz
t0=1/f0;%0.33e-12 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)

semilogx(f/(a0/1e-4),eps);
grid('off')
xlabel('Frequency (Hz)')
ylabel('Relative Permeability \mu_r')

% muinf=1;gamma=.0001;fn=.0001;sigma=-10000;
% f=0:1e-4:(1/3); 
% mur=muinf+(sigma.*fn.*fn)./(-f.*f-1i*gamma*f);
% cf=(1+((f*f0)./(0.2e6)).*((f*f0)./(0.2e6)));
% mur2=mu0+((10000*mu0)./cf)-1i.*cf.*((f*f0)/(0.2e6));
% 
% figure;
% subplot(2,1,1);loglog(f*f0,real(mur));title('Real \mu');
% subplot(2,1,2);loglog(f*f0,imag(mur));title('Imaginary \mu');
% 
% figure;
% subplot(2,1,1);loglog(f*f0,real(mur2)/mu0);title('Real \mu');
% subplot(2,1,2);loglog(f*f0,imag(mur)/mu0);title('Imaginary \mu');
% 
% 
% f=f*f0;mur=mur*mu0;sigma=5e-3;eps=(1*eps0)-1i*(sigma./(2*pi*f));
% gamma=1i.*2.*pi.*f.*sqrt(mur.*eps);
% eta=sqrt(mur./eps);
% Z=gamma.*eta;Y=gamma./eta;
% alpha=real(gamma);beta=imag(gamma);vp=(2.*pi.*f)./beta;
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,abs(Z));title('|Z|');
% subplot(2,1,2);semilogx(f,angle(Z)*(180/pi));title('theta Z');
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,abs(Y));title('|Y|');
% subplot(2,1,2);semilogx(f,angle(Y)*(180/pi));title('theta Y');
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,alpha);title('alpha');
% subplot(2,1,2);semilogx(f,beta);title('beta');
% 
% figure;hold on;
% subplot(2,1,1);semilogx(f,abs(eta));title('\eta');
% subplot(2,1,2);semilogx(f,angle(eta)*(180/pi));title('theta \eta');
% 
% figure;
% semilogx(f,vp);title('vp');
