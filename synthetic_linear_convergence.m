clear all
close all
clc
addpath('tools');

% dimensions of problem
Nt = 1e5;
Ns = 25;
Nb = 6;

% set varibale along each dimension
F = logspace(-3,0,Ns);
B = [5,10,20,25,30,50];  

% set error standard deviation
S = sqrt(0.3);

% set truth
TT = randn(Nt,1);

% loop through experiments
for s = 1:Ns

 % create measurements 
 X = TT + randn(Nt,1)*S;
 Y = TT + randn(Nt,1)*S;
 Z = TT + randn(Nt,1)*S;

 % pull sample
 N = round(Nt*F(s));
 P = randperm(Nt,N);
 X = X(P); Y = Y(P); Z = Z(P); T = TT(P);

 % calculate continuous linear TC stats
 [LE(s,1,1),LE(s,2,1),LE(s,3,1),LI(s,1,1),LI(s,2,1),LI(s,3,1)] = triple_collocation(X,Y,Z);

 % calculate linear continuous truth
 LE(s,1,2) = cov(X-T)/cov(X);
 LE(s,2,2) = cov(Y-T)/cov(Y);
 LE(s,3,2) = cov(Z-T)/cov(Z);
 cc = corrcoef(X,T); LI(s,1,2) = cc(2);
 cc = corrcoef(Y,T); LI(s,2,2) = cc(2);
 cc = corrcoef(Z,T); LI(s,3,2) = cc(2);

 % loop through resolutions
 for b = 1:Nb

  % create bins at resolution
  Bt = linspace(min(T)-1e-6,max(T)+1e-6,B(b));
  Bx = linspace(min(X)-1e-6,max(X)+1e-6,B(b));
  By = linspace(min(Y)-1e-6,max(Y)+1e-6,B(b));
  Bz = linspace(min(Z)-1e-6,max(Z)+1e-6,B(b));

  % calculate nonlinear TC stats
  [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,Bx,By,Bz);

  % bound on total information
  NI(s,1,1,b) = (Ixy+Ixz-Ixyz)/Hx; 
  NI(s,2,1,b) = (Ixy+Iyz-Ixyz)/Hy; 
  NI(s,3,1,b) = (Ixz+Iyz-Ixyz)/Hz; 

  % bound on total error
  NE(s,1,1,b) = 1-(Ixy+Ixz-Ixyz)/Hx;
  NE(s,2,1,b) = 1-(Ixy+Iyz-Ixyz)/Hy;
  NE(s,3,1,b) = 1-(Ixz+Iyz-Ixyz)/Hz;

  % bound on missing information
  NM(s,1,1,b) = (Iyz-Ixyz)/Hx;
  NM(s,2,1,b) = (Ixz-Ixyz)/Hy;
  NM(s,3,1,b) = (Ixy-Ixyz)/Hz;

  % calculate true stats
  [Ixt,Hx,Ht] = mutual_info(X,T,Bx,Bt);
  NI(s,1,2,b) = Ixt/Hx;
  NE(s,1,2,b) = 1-Ixt/Hx;
  NM(s,1,2,b) = (Ht-Ixt)/Hx;

  [Ixt,Hx,Ht] = mutual_info(Y,T,By,Bt);
  NI(s,2,2,b) = Ixt/Hx;
  NE(s,2,2,b) = 1-Ixt/Hx;
  NM(s,2,2,b) = (Ht-Ixt)/Hx;

  [Ixt,Hx,Ht] = mutual_info(Z,T,Bz,Bt);
  NI(s,3,2,b) = Ixt/Hx;
  NE(s,3,2,b) = 1-Ixt/Hx;
  NM(s,3,2,b) = (Ht-Ixt)/Hx;

  % screen report
  [s/Ns,b/Nb]

 end % bin resolution
end % error variance

%% --------- PLOT RESULTS -----------

% grab colors
figure(1); close(1); figure(1);
h = plot(randn(10));
for i = 1:10
 colors(i,:) = h(i).Color;
end
close(1);

binNames(1) = {'linear statistic'};
for b = 1:Nb
 binNames(1+b) = strcat({'bin resolution: '},num2str(round(1/B(b)*1000)/1000),'%');
end

% convergence plots
figure(2); close(2); fig=figure(2);
set(gcf,'color','w','position',[1601,821,875,450]);
%set(gcf,'color','w','position',[1601,821,893,800]);

%subplot(2,1,1)
pdat = squeeze(mean(LE(:,:,1),2));
semilogx(Nt*F,pdat,'-o','color','k','linewidth',2); hold on;
pdat = squeeze(mean(NE(:,:,1,:),2));
semilogx(Nt*F,pdat,'-o'); hold on;
xlabel('number of data','fontsize',16);
ylabel('total error','fontsize',16);
legend(binNames,'location','best');
title('convergence of total error estimates with sample size','fontsize',16);

%subplot(2,1,2)
%pdat = squeeze(mean(LE(:,:,2),2));
%semilogx(Nt*F,pdat,'-o','color','k','linewidth',2); hold on;
%pdat = squeeze(mean(NE(:,:,2,:),2));
%semilogx(Nt*F,pdat,'-o'); hold on;
%xlabel('number of data','fontsize',16);
%ylabel('total error','fontsize',16);
%legend(binNames,'location','best');
%title('convergence of total error (true) with sample size','fontsize',16);

fname = 'figures/Figure5_LinearConvergence';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

