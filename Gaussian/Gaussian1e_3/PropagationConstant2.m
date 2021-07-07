clc;clear all;
%hold on;
%%
f=fopen('FieldEvolutionIn.txt');
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
    
    data{in}=sscanf(text{in},'%f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f');
    dis=data{in}(2)+13;
%    dis=round(data{in}(2)*10+61);
    Hx(data{in}(1)-Prev,dis)=data{in}(3)+1i*data{in}(4);
    Hy(data{in}(1)-Prev,dis)=data{in}(5)+1i*data{in}(6);
    Hz(data{in}(1)-Prev,dis)=data{in}(7)+1i*data{in}(8);
    
    Bx(data{in}(1)-Prev,dis)=data{in}(9)+1i*data{in}(10);
    By(data{in}(1)-Prev,dis)=data{in}(11)+1i*data{in}(12);
    Bz(data{in}(1)-Prev,dis)=data{in}(13)+1i*data{in}(14);
    
    Ex(data{in}(1)-Prev,dis)=data{in}(15)+1i*data{in}(16);
    Ey(data{in}(1)-Prev,dis)=data{in}(17)+1i*data{in}(18);
    Ez(data{in}(1)-Prev,dis)=data{in}(19)+1i*data{in}(20);
   
    Dx(data{in}(1)-Prev,dis)=data{in}(21)+1i*data{in}(22);
    Dy(data{in}(1)-Prev,dis)=data{in}(23)+1i*data{in}(24);
    Dz(data{in}(1)-Prev,dis)=data{in}(25)+1i*data{in}(26);
    
    l=fgetl(f);
    in=in+1;
end

epsr=1;
a0=1e-2;%0.1mm
c0=3e8;%2.99792458e8;%Speed of Light (m/s)
f0=c0/a0;%3000GHz
t0=1/f0;%0.33e-12 (s)
mu0=4*pi*(1e-7);% (H/m)
eps0=8.854187817e-12;% (F/m)
I0=1; %(A)
E0=I0/(a0*eps0*c0);%Electric Field
D0=I0/(a0*c0);%Electric Displacement Field
B0=I0/(a0*eps0*c0*c0);%Magnetic Field
H0=I0/(a0);%Magnetizing Field
fmin=1e10;
fmax=5e10;

Lmax=1024;
L2=Lmax;

zi=14;
zo=15;
subplot(2,2,1)
plot(((1:Lmax)*t0)/50,real(Hx(1:Lmax,zi))*H0);hold on;
plot(((1:Lmax)*t0)/50,real(Hy(1:Lmax,zi)*H0),'r');
xlabel('Time(s)');ylabel('Magnitude (A/m)');title('Evolution of Fields at Source End');
legend('Hx','Hy');
%axis([0 Lmax*t0 -1000 1000]);
subplot(2,2,2)
hold on
for tind=1:Lmax
plot([real(Hx(tind,zi)) real(Hx(tind+1,zi))]*H0,[real(Hy(tind,zi)) real(Hy(tind+1,zi))]*H0)
%axis([-0.1 0.1 -0.1 0.1]);
end
xlabel('Hx (A/m)');ylabel('Hy (A/m)');title('Evolution of Fields at Source End');axis('equal');
%axis([-1000 1000 -1000 1000]);
subplot(2,2,3)
plot(((1:Lmax)*t0)/50,real(Hx(1:Lmax,zo))*H0);hold on;
plot(((1:Lmax)*t0)/50,real(Hy(1:Lmax,zo)*H0),'r');
xlabel('Time(s)');ylabel('Magnitude (A/m)');
%axis([0 Lmax*t0 -200 200]);
legend('Hx','Hy');
subplot(2,2,4)
hold on
for tind=1:Lmax
plot([real(Hx(tind,zo)) real(Hx(tind+1,zo))]*H0,[real(Hy(tind,zo)) real(Hy(tind+1,zo))]*H0)
%axis([-0.1 0.1 -0.1 0.1]);
end
xlabel('Hx (A/m)');ylabel('Hy (A/m)');axis('equal');
%axis([-200 200 -200 200]);


T=t0/50;
Fs=1/T;
L=Lmax;

dobs=5+1;
obs=14;
L=2^nextpow2(L);
f=Fs/2*linspace(0,1,L/2+1);
Heffi=eff(Hx(1:L,obs),Hy(1:L,obs));
Heffo=eff(Hx(1:L,obs+dobs),Hy(1:L,obs+dobs));
FHxi=(fft(Heffi,L));
FHxo=(fft(Heffo,L));
Gammap=(log(FHxo./FHxi))/(-(dobs/1)*(a0));


dobs=-5-1;
obs=12;
L=2^nextpow2(L);
f=Fs/2*linspace(0,1,L/2+1);
Heffi=eff(Hx(1:L,obs),Hy(1:L,obs));
Heffo=eff(Hx(1:L,obs+dobs),Hy(1:L,obs+dobs));
FHxi=(fft(Heffi,L));
FHxo=(fft(Heffo,L));
Gamman=(log(FHxo./FHxi))/((dobs/1)*(a0));


figure;
subplot(2,1,1);hold on;
plot(f(1:L/2+1),real(Gammap(1:L/2+1)));hold on;
plot(f(1:L/2+1),real(Gamman(1:L/2+1)),'r');
ylabel('\alpha (Np.m^-^1)');xlim([fmin fmax]);
subplot(2,1,2)
plot(f(1:L/2+1),imag(Gammap(1:L/2+1)));hold on;
plot(f(1:L/2+1),imag(Gamman(1:L/2+1)),'r');
ylabel('\beta (rad.m^-^1)');xlabel('Frequency (Hz)');xlim([fmin fmax]);


obs=14;
Eeffo=eff(Ex(1:L,obs),Ey(1:L,obs));
Heffo=eff(Hx(1:L,obs),Hy(1:L,obs));
NFFT=2^nextpow2(L);
f=Fs/2*linspace(0,1,NFFT/2+1);
FEx=fft(Eeffo,NFFT);
FHx=fft(Heffo,NFFT);
Zp=FEx./FHx;
Zp=Zp*377;

obs=12;
Eeffo=eff(Ex(1:L,obs),Ey(1:L,obs));
Heffo=eff(Hx(1:L,obs),Hy(1:L,obs));
NFFT=2^nextpow2(L);
f=Fs/2*linspace(0,1,NFFT/2+1);
FEx=fft(Eeffo,NFFT);
FHx=fft(Heffo,NFFT);
Zn=FEx./FHx;
Zn=Zn*377;


figure;
subplot(2,1,1)
plot(f(1:NFFT/2+1),abs(Zp(1:NFFT/2+1)));hold on;
plot(f(1:NFFT/2+1),abs(Zn(1:NFFT/2+1)),'r');
ylabel('|\eta| (Ohm)');xlim([fmin fmax]);
hold on;
subplot(2,1,2)
semilogx(f(1:NFFT/2+1),abs(angle(Zp(1:NFFT/2+1))*(180/pi)));hold on;
semilogx(f(1:NFFT/2+1),abs(angle(Zn(1:NFFT/2+1))*(180/pi)),'r');
ylabel('\Theta \eta (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);


% Zimag=abs(imag(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
% Zreal=abs(real(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zabsp=abs((Gammap(1:NFFT/2+1).*Zp(1:NFFT/2+1)));
Zanglep=abs(angle((Gammap(1:NFFT/2+1).*Zp(1:NFFT/2+1)))*(180/pi));
% Gm=abs(real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
% Rm=abs((Gm')*(-1).*(2*pi*f));%Reluctance
% XCm=abs(imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Yabsp=abs((Gammap(1:NFFT/2+1)./Zp(1:NFFT/2+1)));
Yanglep=(angle((Gammap(1:NFFT/2+1)./Zp(1:NFFT/2+1)))*(180/pi));

% Zimag=abs(imag(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
% Zreal=abs(real(Gamma(1:NFFT/2+1).*Z(1:NFFT/2+1)));
Zabsn=abs((Gamman(1:NFFT/2+1).*Zn(1:NFFT/2+1)));
Zanglen=abs(angle((Gamman(1:NFFT/2+1).*Zn(1:NFFT/2+1)))*(180/pi));
% Gm=abs(real(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
% Rm=abs((Gm')*(-1).*(2*pi*f));%Reluctance
% XCm=abs(imag(Gamma(1:NFFT/2+1)./Z(1:NFFT/2+1)));
Yabsn=abs((Gamman(1:NFFT/2+1)./Zn(1:NFFT/2+1)));
Yanglen=(angle((Gamman(1:NFFT/2+1)./Zn(1:NFFT/2+1)))*(180/pi));

figure;
subplot(2,1,1)
plot(f,(Zabsp(1:NFFT/2+1)));hold on;
plot(f,(Zabsn(1:NFFT/2+1)),'r');
ylabel('|Z| (Ohm/m)');xlim([fmin fmax]);title('Transverse Impedance Z');
hold on;

subplot(2,1,2);
semilogx(f,abs(Zanglep(1:NFFT/2+1)));hold on;
semilogx(f,abs(Zanglen(1:NFFT/2+1)),'r');
ylabel('\theta Z (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
hold on;

% subplot(3,1,1);
% plot(f(1:NFFT/2+1),(Gm(1:NFFT/2+1)));
% ylabel('Conductance GL (S/m)');
% xlabel('Frequency (Hz)');
% % axis([0 5e11 -8.5 -7.5])
% 
% subplot(3,1,2);
% plot(f(1:NFFT/2+1),(Rm(1:NFFT/2+1)));
% ylabel('Reluctance Rmskin (1/H.m)');
% xlabel('Frequency (Hz)');
% % axis([0 5e11 0 3e13])
% 
% subplot(3,1,3);
% plot(f(1:NFFT/2+1),(XCm(1:NFFT/2+1)));
% ylabel('Susceptance XCL (S/m)');
% xlabel('Frequency (Hz)');
% % axis([0 5e11 0 1])

figure;
subplot(2,1,1)
plot(f,(Yabsp(1:NFFT/2+1)));hold on;
plot(f,(Yabsn(1:NFFT/2+1)),'r');
ylabel('Y (1/Ohm.m)');xlim([fmin fmax]);

subplot(2,1,2);
semilogx(f,abs(Yanglep(1:NFFT/2+1)));hold on;
semilogx(f,abs(Yanglen(1:NFFT/2+1)),'r');
ylabel('\theta Y (Degree)');xlabel('Frequency (Hz)');xlim([fmin fmax]);
