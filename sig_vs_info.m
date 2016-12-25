clear all
close all
clc

addpath('tools');

N = 1e5;
M = 50;

sig = linspace(0,sqrt(5),M);

B = -50:0.1:50;

for m = 1:M
 T = randn(N,1);
 X = T+randn(N,1)*sig(m);
 [Ixt(m),Hx(m),Ht(m)] = mutual_info(X,T,B,B); 
 cc = corrcoef(X,T); rho(m) = cc(2);
 m/M
end

H = (Hx-Ixt)./Hx;
I = Ixt./Hx;
S = sig.^2;

figure(1); close(1); fig = figure(1);
set(gcf,'color','w','position',[3200,450,750,1000]);

subplot(2,1,1)
[ax,h1,h2] = plotyy(I,S,I,rho); 
ylabel(ax(1),'error variance','fontsize',18);
set(ax(1),'YColor','k');
h1.LineWidth = 4;
h1.Color = [0,0,0];
ylabel(ax(2),'linear correlation','fontsize',18);
set(ax(2),'YColor','k');
h2.LineWidth = 4;
h2.Color = [0.4,0.4,0.4];
xlabel('information ratio','fontsize',18);
title('Error Variance vs. Measurement Information','fontsize',16)

subplot(2,1,2)
[ax,h1,h2] = plotyy(H,S,H,rho); 
ylabel(ax(1),'error variance','fontsize',18);
set(ax(1),'YColor','k');
h1.LineWidth = 4;
h1.Color = [0,0,0];
ylabel(ax(2),'linear correlation','fontsize',18);
set(ax(2),'YColor','k');
h2.LineWidth = 4;
h2.Color = [0.4,0.4,0.4];
xlabel('total error fraction','fontsize',18);
title('Error Variance vs. Measurement Error','fontsize',16)

fname = 'figures/Figure2_sigVinfo';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);



[sig',Ixt'./Hx',rho']


