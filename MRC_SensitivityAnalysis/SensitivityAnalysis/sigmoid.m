p = [0:150];
HM = 3;
Hm = .5;
k = 8;
p2 = 90;
x1 = (HM-Hm)*p2.^k./(p.^k+p2^k)+Hm;
x2 = (HM-Hm)*p.^k./(p.^k+p2^k)+Hm;
figure(1);
h=plot(p,x1,p,x2)
set(h,'linewidth',4);

figure(1);
print -depsc2 sigmoids.eps
