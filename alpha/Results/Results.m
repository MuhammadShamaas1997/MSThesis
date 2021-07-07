Alpha=[1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2];
YL=[0.009902 0.009872 0.009841 0.009872 0.01065 0.0583 0.46];
ZT=[3.241e5 3.244e5 3.244e5 3.159e5 2.841e5 6.678e4 1.393e4];
ThetaZT=[16.89 17.52 18.31 25.04 33.56 67.44 74.69];
Eta=[5271 5732 5742 5657 5163 1070 174];
ThetaEta=[64.93 65.62 66.51 73.93 83.39 121.6 90.52];
alpha=[37.88 37.79 37.67 36.72 35.49 36.56 77.01];
beta=[-42.12 -42.12 -42.12 -42.08 -42.04 -50.56 -21.83];

semilogx(Alpha,YL,'-o')
figure;semilogx(Alpha,ZT,'-o')
figure;semilogx(Alpha,ThetaZT,'-o')
figure;semilogx(Alpha,Eta,'-o')
figure;semilogx(Alpha,ThetaEta,'-o')
figure;semilogx(Alpha,alpha,'-o')
figure;semilogx(Alpha,beta,'-o')
