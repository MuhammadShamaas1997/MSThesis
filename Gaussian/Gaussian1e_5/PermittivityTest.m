clc;clear all;

fi=fopen('Permeability.txt');
l=fgetl(fi);
in=1;
while ischar(l)
%for kj=1:81
    %%disp(l);
    text{in}=l;
    data{in}=sscanf(text{in},'%f (%f,%f)');
    f(in)=(data{in}(1));
    eps(in)=data{in}(2)+1i*data{in}(3);
    
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
I0=1;
E0=I0/(a0*eps0*c0);
D0=I0/(a0*c0);
B0=I0/(a0*eps0*c0*c0);
H0=I0/(a0);
sigmaD0=(epsr*eps0*c0)/a0;
J0=I0/(a0*a0);
u0=(I0*I0)/(eps0*c0*c0*a0*a0);
S0=(I0*I0)/(eps0*c0*a0*a0);
Sc0=1/(c0);
sig0=-20;


figure;
fdata=[1000,2000,3000,4000,8000,10000,12000,15000,20000,25000,30000,40000,50000,60000,75000,80000,96000,100000,150000,200000,300000,400000,500000,600000,800000,960000,1000000,1137137,1274275,1438450,1623777,1832981,2084351,2335721,2636651,2976351,3359818,3792690,4281332,4868464,5455595,6158482,6951928,7847600,8923800,10000000];
muredata=[10645.6132900000,10615.7137400000,10600.4558900000,10588.2058700000,10566.9610700000,10562.5965000000,10550.3503700000,10534.2909700000,10511.8581200000,10483.3333400000,10455.5688100000,10381.6845200000,10272.9184600000,10106.1882500000,9714.59266900000,9554.00022600000,8914.62151900000,8809.17242500000,6859.02588800000,5356.50053100000,3466.60912500000,2509.42420600000,1981.99317800000,1653.79743600000,1259.89804500000,1061.94435800000,1039.71574800000,921.731119300000,830.272918200000,736.042665000000,651.396975700000,577.613199300000,501.648967800000,445.436779500000,389.427157200000,340.949677100000,297.778927600000,260.157345700000,226.241196700000,196.013417100000,172.718377000000,149.043074900000,129.745046100000,112.942358600000,97.9382173800000,86.5638240400000];
muimdata=[253.852268500000,226.074739600000,240.928212700000,269.680112600000,420.349094000000,501.445173000000,580.251521400000,703.380303600000,910.425579500000,1120.16008100000,1332.71255700000,1775.31979800000,2241.36105700000,2719.66949700000,3410.37939800000,3618.92417300000,4304.31495800000,4315.37839400000,5129.34009300000,5172.35110400000,4575.54022900000,3898.79388100000,3374.31137400000,2983.65236200000,2453.70814500000,2168.65609600000,2103.25394200000,1931.01903200000,1787.88868700000,1647.75656600000,1514.22011300000,1390.88856500000,1270.26340100000,1169.77709400000,1067.60990700000,974.786626800000,888.046742400000,806.937314600000,733.896914500000,660.333778100000,601.248498000000,543.010068900000,490.624234300000,442.324718400000,394.992484200000,356.844963000000];

subplot(2,1,1);hold on;
semilogx(f/(a0/(1e-2)),real(eps));
semilogx(fdata,muredata,'r');xlim([1e4 .2e7]);
title('Relative Permeability \mu_r')
ylabel('Real \mu_r')

subplot(2,1,2);hold on;
semilogx(f/(a0/(1e-2)),imag(eps));
semilogx(fdata,-muimdata,'r');xlim([1e4 .2e7]);

xlabel('Frequency (Hz)')
ylabel('Imaginary')

figure;
semilogx(f/(a0/(1e-2)),real(eps));
hold on;
semilogx(fdata,muredata,'r');
xlim([1e4 1e7]);
title('Relative Permeability \mu_r')
semilogx(f/(a0/(1e-2)),-imag(eps),'-.');
semilogx(fdata,muimdata,'-.');
legend('Real (Simulation)','Real (Datasheet)','Imaginary (Simulation)','Imaginary (Datasheet)');
% muinf=1;gamma=.01/(8*4*pi*pi*pi*pi*pi);fn=0.01/(4*pi*pi);
% f=0:1e-7:(1/3);
% fi=f*(4*pi*pi);
% sigma=-100*(4*pi*pi);
% kwi=(fn.*fn-fi.*fi-1i.*gamma.*fi);
% b0=1.0;
% kwi2=(fi.*fi.*b0.*b0)./kwi;
% mur=muinf+(sigma.*fn.*fn)./(kwi-kwi2);
b0=0;
fn=(5.8)*(1e-5);
sigma=1.0562e4;
muinf=1;gamma=-(.33e-2)/(2*pi);
f=0:1e-7:(1/3);
fi=f;
kwi=(fn.*fn-fi.*fi-1i.*gamma.*fi);
kwi2=(fi.*fi.*b0.*b0)./kwi;
mur=(sigma.*fn.*fn.*kwi)./(kwi.*kwi-fi.*fi.*b0.*b0);

figure;
subplot(2,1,1);semilogx(f*f0,real(mur));title('Formula');
hold on;semilogx(fdata,muredata,'r');
xlim([1e4 1e7]);
subplot(2,1,2);semilogx(f*f0,imag(mur));xlim([10e5 1e9]);
hold on;semilogx(fdata,-muimdata,'r');
xlim([1e4 1e7]);

f=1e3:1e3:1e9;
nomf=((f)./(0.2e6));
cf=(1+nomf.*nomf);
mur2=mu0+((10000*mu0)./cf)-(1i.*nomf.*(10000*mu0))./cf;


% figure;
% subplot(2,1,1);hold on;
% f=0:1e-7:(1/30);
% loglog(f*f0,real(mur));title('Formula');xlim([1e3 100e5]);
% f=1e3:1e3:1e9;
% loglog(f,real(mur2)/mu0,'r');title('Literature');xlim([1e3 100e5]);
% subplot(2,1,2);hold on;
% f=0:1e-7:(1/30);
% loglog(f*f0,imag(mur));xlim([1e3 100e5]);
% f=1e3:1e3:1e9;
% loglog(f,imag(mur2)/mu0,'r');xlim([1e3 100e5]);

b0=0;
fn=(5.8)*(1e-5);
sigma=1.0562e4;
muinf=1;gamma=-(.33e-2)/(2*pi);
f=0:1e-7:(1/3);
fi=f;
kwi=(fn.*fn-fi.*fi-1i.*gamma.*fi);
kwi2=(fi.*fi.*b0.*b0)./kwi;
mur=(sigma.*fn.*fn.*kwi)./(kwi.*kwi-fi.*fi.*b0.*b0);


alpha=0.00001;
fn=1;
sigma=.1;%1.0562e4;
muinf=1;gamma=(.001)/(2*pi);
f=0.99:1e-4:1.01;
fi=f;
kwi=(fn-1i.*fi.*alpha);
kwi2=fi+1i*gamma;
mur=(sigma.*kwi)./(kwi.*kwi-kwi2.*kwi2);

figure;
subplot(2,1,1);
hold on;
plot(f*f0,abs(mur));ylabel('\chi_1_1');xlabel('Frequency')
title('Gyromagnetic Susceptibility');
plot(f*f0,imag(mur),'r');legend('Real','Imaginary');
%hold on;semilogx(fdata,-muimdata,'r');
%xlim([1e4 1e7]);



% alpha=0.00001;
% fn=0.001;
% sigma=.1;%1.0562e4;
% muinf=1;gamma=(.0001)/(2*pi);
% f=0:1e-7:0.002;
fi=f;
kwi=(fn-1i.*fi.*alpha);
kwi2=fi+1i*gamma;
mur=(sigma.*kwi2)./(kwi.*kwi-kwi2.*kwi2);

subplot(2,1,2);
hold on;
plot(f*f0,real(mur));ylabel('\chi_1_2');xlabel('Frequency')
%title('Gyromagnetic Susceptibility');
plot(f*f0,imag(mur),'r');legend('Real','Imaginary');
%hold on;semilogx(fdata,-muimdata,'r');
%xlim([1e4 1e7]);
