clear all
close all
clc

% dimensions of problem
Nt = 1e6;
Ns = 7;
Nb = 1;

% set varibale along each dimension
S = linspace(1e-3,sqrt(0.5),Ns);
B = [25];%linspace(5,100,Nb);

% set truth
T = rand(Nt,1);

% loop through experiments
for s = 1:Ns
 
 % create measurements 
% X = T  + randn(Nt,1)*S(s);
% Y = T .* rand(Nt,1)*S(s);
% Z = T .* (1+randn(Nt,1)*S(s));
 X = T  + randn(Nt,1)*S(s);
 Y = T .* randn(Nt,1)*S(s);
 Z = T .* (1+randn(Nt,1)*S(s));

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
%  By = logspace(log10(1),log10(max(Y))+1e-6,round(B(b)/2)); By = [-fliplr(By),By];
  Bz = linspace(min(Z)-1e-6,max(Z)+1e-6,B(b));

  % calculate nonlinear TC stats
  [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,Bx,By,Bz);

  % bound on total info mation
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

% scatter plots
figure(1); close(1); fig=figure(1);
set(gcf,'color','w','position',[3200,450,725,1200]);
%set(gcf,'color','w','position',[3200,450,450,450]);

subplot(2,1,1)
plot(squeeze(LE(:,:,2)),squeeze(LE(:,:,1)),'-o','linewidth',2); hold on;
plot(squeeze(LI(:,:,2)),squeeze(LI(:,:,1)),'--+','linewidth',2); hold on;
plot([0,1],[0,1],'k--')
grid on; 
title('Linear Triple Collocation','fontsize',16);
xlabel('true statistic','fontsize',18);
ylabel('estimated statistic','fontsize',18);
legend('total error: additive noise','total error: multiplicative noise','total error: heteroscedastic noise',...
       'correlation: additive noise','correlation: multiplicative noise','correlation: heteroscedastic noise','location','nw');
%axis([-0.5,2,-0.5,2]);
axis([0,1,0,1]);

subplot(2,1,2)
plot(squeeze(NE(:,:,2,ceil(Nb/3))),squeeze(NE(:,:,1,ceil(Nb/3))),'-o','linewidth',2); hold on;
plot(squeeze(NI(:,:,2,ceil(Nb/3))),squeeze(NI(:,:,1,ceil(Nb/3))),'--+','linewidth',2); hold on;
%plot(squeeze(NM(:,:,2,ceil(Nb/3))),squeeze(NM(:,:,1,ceil(Nb/3))),':s','linewidth',2); hold on;
plot([0,1],[0,1],'k--')
grid on; 
title('Nonlinear Triple Collocation','fontsize',16);
xlabel('true statistic','fontsize',18);
ylabel('estimated statistic','fontsize',18);
legend('total error: additive noise','total error: multiplicative noise','total error: heteroscedastic noise',...
       'total info : additive noise','total info : multiplicative noise','total info : heteroscedastic noise','location','nw');
axis([0,1,0,1]);

fname = 'figures/Figure6_NonlinearSynthetic';
img = getframe(gcf);
imwrite(img.cdata, [fname,'.png']);











