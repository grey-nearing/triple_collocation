clear all
close all
clc
addpath('tools');

% dimensions 
Ns = 100;
Nb = 3;

% read data
data = importdata('data/wgew_precip/gauge_data.txt');
date = data.data(:,1:3);
data = data.data(:,4:end);

% number of samples and gauges
[Nt,Ng] = size(data);

% truth
T = mean(data,2);

% set varibale along each dimension
B = [50,20,10];%linspace(5,1000,Nb);
plotB = 2;

% random permutations
perms = zeros(Ns,3);
%for m = 1:Ns
% perms(m,:) = randperm(Ng,3);
%end

% loop through experiments
for s = 1:Ns

 indep = ones(3,1);
 while max(indep)>0.25
  
  % get a sample of gauges
  perms(s,:) = randperm(Ng,3);
 
  % create measurements 
  X = data(:,perms(s,1));
  Y = data(:,perms(s,2));
  Z = data(:,perms(s,3));

  Bt = linspace(min(T)-1e-6,max(T)+1e-6,B(3));
  Bx = linspace(min(X)-1e-6,max(X)+1e-6,B(3));
  By = linspace(min(Y)-1e-6,max(Y)+1e-6,B(3));
  Bz = linspace(min(Z)-1e-6,max(Z)+1e-6,B(3));

  % measure independence
  [Ixyt,Ixy,Ixt,Iyt,Hx,Hy,Ht] = mutual_info_3(X,Y,T,Bx,By,Bt);
  [Ixtz,Ixt,Ixz,Itz,Hx,Ht,Hz] = mutual_info_3(X,T,Z,Bx,Bt,Bz);
  [Ityz,Ity,Itz,Iyz,Ht,Hy,Hz] = mutual_info_3(T,Y,Z,Bt,By,Bz);

  indep(1) = (Ixy-Ixyt)/Ixy;
  indep(2) = (Ixz-Ixtz)/Ixz;
  indep(3) = (Iyz-Ityz)/Iyz;

  [s,indep',Ixyt,Ixtz,Ityz]

 end

% % remove grandmas
% I = find(any([X,Y,Z]'<=0));
% X(I) = []; Y(I) = []; Z(I) = []; T(I) = [];

 % log transform
 lX = max(-6,log(X));
 lY = max(-6,log(Y));
 lZ = max(-6,log(Z));
 lT = max(-6,log(T));

 % calculate continuous linear TC stats
 [LE(s,1,1),LE(s,2,1),LE(s,3,1),LI(s,1,1),LI(s,2,1),LI(s,3,1)] = triple_collocation(X,Y,Z);
 [GE(s,1,1),GE(s,2,1),GE(s,3,1),GI(s,1,1),GI(s,2,1),GI(s,3,1)] = triple_collocation(lX,lY,lZ);

 % calculate linear continuous truth
 LE(s,1,2) = cov(X-T)/cov(X);
 LE(s,2,2) = cov(Y-T)/cov(Y);
 LE(s,3,2) = cov(Z-T)/cov(Z);
 cc = corrcoef(X,T); LI(s,1,2) = cc(2);
 cc = corrcoef(Y,T); LI(s,2,2) = cc(2);
 cc = corrcoef(Z,T); LI(s,3,2) = cc(2);

 % linear log-transformed truth
 GE(s,1,2) = cov(lX-lT)/cov(lX);
 GE(s,2,2) = cov(lY-lT)/cov(lY);
 GE(s,3,2) = cov(lZ-lT)/cov(lZ);
 cc = corrcoef(lX,lT); GI(s,1,2) = cc(2);
 cc = corrcoef(lY,lT); GI(s,2,2) = cc(2);
 cc = corrcoef(lZ,lT); GI(s,3,2) = cc(2);

 % loop through resolutions
 for b = 1:Nb

  % create bins at resolution
%  Bt = [0,logspace(log10(0.001),log10(max(T)),B(b))]; Bt(end) = Bt(end)+1e-6;
%  Bx = [0,logspace(log10(0.001),log10(max(X)),B(b))]; Bx(end) = Bx(end)+1e-6;
%  By = [0,logspace(log10(0.001),log10(max(Y)),B(b))]; By(end) = By(end)+1e-6;
%  Bz = [0,logspace(log10(0.001),log10(max(Z)),B(b))]; Bz(end) = Bz(end)+1e-6;
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

  % observable correlations
  CC(s,1,2,b) = Ixy/Hx; CC(s,1,3,b) = Ixz/Hx;
  CC(s,2,1,b) = Ixy/Hy; CC(s,2,3,b) = Iyz/Hy;
  CC(s,3,1,b) = Ixz/Hz; CC(s,3,2,b) = Iyz/Hz;

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

  % measure independence
  [Ixyt,Ixy,Ixt,Iyt,Hx,Hy,Ht] = mutual_info_3(X,Y,T,Bx,By,Bt);
  [Ixtz,Ixt,Ixz,Itz,Hx,Ht,Hz] = mutual_info_3(X,T,Z,Bx,Bt,Bz);
  [Ityz,Ity,Itz,Iyz,Ht,Hy,Hz] = mutual_info_3(T,Y,Z,Bt,By,Bz);
  IND(s,1,b) = (Ixy-Ixyt)/Ixy;
  IND(s,2,b) = (Ixz-Ixtz)/Ixz;
  IND(s,3,b) = (Iyz-Ityz)/Iyz;

  % screen report
  [s/Ns,b/Nb]

 end % bin resolution
end % error variance

save('results/wgew_raw_stats.mat');
%load('wgew_raw_stats.mat');

% linear stats
let = LE(:,:,2); let = let(:);
lee = LE(:,:,1); lee = lee(:);
lit = LI(:,:,2); lit = lit(:);
lie = LI(:,:,1); lie = lie(:);
let = LE(:,:,2); let = let(:);
lee = LE(:,:,1); lee = lee(:);
lit = LI(:,:,2); lit = lit(:);
lie = LI(:,:,1); lie = lie(:);
get = GE(:,:,2); get = get(:);
gee = GE(:,:,1); gee = gee(:);
git = GI(:,:,2); git = git(:);
gie = GI(:,:,1); gie = gie(:);

% loss functions (bias)
le_bias = mean(abs(let-lee));
li_bias = mean(abs(lit-lie));
ge_bias = mean(abs(get-gee));
gi_bias = mean(abs(git-gie));

% loss functions (var)
le_var = mean((let-lee).^2);
li_var = mean((lit-lie).^2);
ge_var = mean((get-gee).^2);
gi_var = mean((git-gie).^2);

% loss functions (mse) 
le_mse = le_bias^2+le_var;
li_mse = li_bias^2+li_var;
ge_mse = ge_bias^2+ge_var;
gi_mse = gi_bias^2+gi_var;

% loss funciton correlation
cc = corrcoef(let,lee); le_cor = cc(2);
cc = corrcoef(lit,lie); li_cor = cc(2);
cc = corrcoef(get,gee); ge_cor = cc(2);
cc = corrcoef(git,gie); gi_cor = cc(2);

% nonlinsear stats
for b = 1:Nb

 net = NE(:,:,2,b); net = net(:);
 nee = NE(:,:,1,b); nee = nee(:);
 nit = NI(:,:,2,b); nit = nit(:);
 nie = NI(:,:,1,b); nie = nie(:);

 ne_bias(b) = mean(abs(net-nee));
 ni_bias(b) = mean(abs(nit-nie));

 % loss functions (variance)
 ne_var(b) = mean((net-nee).^2);
 ni_var(b) = mean((nit-nie).^2);

 % loss functions (mse)
 ne_mse(b) = ne_bias(b)^2+ne_var(b);
 ni_mse(b) = ni_bias(b)^2+ni_var(b);

 % correlation
 cc = corrcoef(net,nee); ne_cor(b) = cc(2);
 cc = corrcoef(nit,nie); ni_cor(b) = cc(2);

 % risk significance tests
 [h,p] = ttest(abs(net-nee),abs(let-lee),'Alpha',1e-3);
 [b,mean(abs(net-nee)),mean(abs(let-lee)),h,p]

 [h,p] = ttest(abs(nit-nie),abs(lit-lie),'Alpha',1e-3);
 [b,mean(abs(nit-nie)),mean(abs(lit-lie)),h,p]

 [h,p] = ttest(abs(net-nee),abs(get-gee),'Alpha',1e-3);
 [b,mean(abs(net-nee)),mean(abs(get-gee)),h,p]

 [h,p] = ttest(abs(nit-nie),abs(git-gie),'Alpha',1e-3);
 [b,mean(abs(nit-nie)),mean(abs(git-gie)),h,p]

end

% screen report
error_bias = [ne_bias,le_bias,ge_bias]
info_bias  = [ni_bias,li_bias,gi_bias]
error_var  = [ne_var ,le_var ,ge_var ]
info_var   = [ni_var ,li_var ,gi_var ]
error_mse  = [ne_mse ,le_mse ,ge_mse ]
info_mse   = [ni_mse ,li_mse ,gi_mse ]
error_cor  = [ne_cor ,le_cor ,ge_cor ]
info_cor   = [ni_cor ,li_cor ,gi_cor ]

%grab colors for ploting
figure(100); h = plot(randn(10,10));

%%% independence plots
%figure(10); close(10); figure(10);
%set(gcf,'color','w','position',[3200,450,450,450]);
%
%pdat = squeeze(IND);
%hist(pdat(:));
%set(get(gca,'child'),'facecolor','k','edgecolor','w');
%xlabel('info frac [~]','fontsize',18);
%ylabel('frequency','fontsize',18);
%title('Consitional Sample Dependence','fontsize',16);

% scatter plots
figure(1); close(1); figure(1);
set(gcf,'color','w','position',[3200,450,450,1200]);

subplot(3,1,1)
h1 = plot(squeeze(LE(:,:,2)),squeeze(LE(:,:,1)),'o','color',h(4).Color); hold on;
h2 = plot(squeeze(LI(:,:,2)),squeeze(LI(:,:,1)),'+','color',h(2).Color); hold on;
plot([0,1],[0,1],'k--')
grid on; 
title('Linear TC','fontsize',16);
xlabel('true statistic','fontsize',18);
ylabel('estimated statistic','fontsize',18);
legend([h1(1),h2(1)],'total error','correlation','location','nw');
axis([0,1,0,1]);

subplot(3,1,2)
h1 = plot(squeeze(GE(:,:,2)),squeeze(GE(:,:,1)),'o','color',h(4).Color); hold on;
h2 = plot(squeeze(GI(:,:,2)),squeeze(GI(:,:,1)),'+','color',h(2).Color); hold on;
plot([0,1],[0,1],'k--')
grid on; 
title('Linear Log-Transformed TC','fontsize',16);
xlabel('true statistic','fontsize',18);
ylabel('estimated statistic','fontsize',18);
legend([h1(1),h2(1)],'total error','correlation','location','nw');
axis([0,1,0,1]);

subplot(3,1,3)
h1 = plot(squeeze(NE(:,:,2,plotB)),squeeze(NE(:,:,1,plotB)),'o','color',h(4).Color); hold on;
h2 = plot(squeeze(NI(:,:,2,plotB)),squeeze(NI(:,:,1,plotB)),'+','color',h(2).Color); hold on;
plot([0,1],[0,1],'k--')
grid on; 
title('Nonlinear TC','fontsize',16);
xlabel('true statistic','fontsize',18);
ylabel('estimated statistic','fontsize',18);
legend([h1(1),h2(1)],'total error','total information','location','nw');
axis([0,1,0,1]);

fname = 'figures/Figure5_wgewGaugeOnly';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

asdf


% sensitivity plots
for b = 1:Nb
 ii = NI(:,:,2,b); ii = ii(:); TI(:,b) = ii;
 ii = NI(:,:,1,b); ii = ii(:); EI(:,b) = ii;
 ee = NE(:,:,2,b); ee = ee(:); TE(:,b) = ee;
 ee = NE(:,:,1,b); ee = ee(:); EE(:,b) = ee;
 cc = CC(:,:,:,b); cc = cc(:); CI(:,b) = cc;
end

%for b = 1:Nb
% ii = (NI(:,:,2,b)-NI(:,:,1,b))./NI(:,:,2,b); ii = ii(:); TI(:,b) = ii;
% ee = (NE(:,:,2,b)-NE(:,:,1,b))./NE(:,:,2,b); ee = ee(:); TE(:,b) = ee;
%end

figure(2); close(2); figure(2);
set(gcf,'color','w','position',[3200,450,725,725]);

%h1 = errorbar(B,mean(TI),std(TE),'o-','linewidth',2,'color',h(4).Color); hold on;
%h2 = errorbar(B,mean(TE),std(TE),'o--','linewidth',2,'color',h(2).Color); hold on;
%h1 = errorbar(B,mean(TI),std(TI),'o-','linewidth',2,'color',h(4).Color); hold on;
%h2 = errorbar(B,mean(EE),std(EE),'o--','linewidth',2,'color',h(1).Color); hold on;
%h3 = errorbar(B,mean(TI),std(TI),'o-','linewidth',2,'color',h(2).Color); hold on;
%h4 = errorbar(B,mean(EI),std(EI),'o--','linewidth',2,'color',h(3).Color); hold on;
h5 = plot(B,CI,'o-','linewidth',2,'color',h(2).Color); hold on;
h1 = errorbar(B,mean(TI),std(TI),'o-','linewidth',2,'color',h(4).Color); hold on;


%legend([h1,h2,h3,h4],'true error','estimated error','true info','estimated info','location','nw');
legend([h1,h5(1)],'true total info','observed cross-info','location','nw');
xlabel('number of bins','fontsize',16);
ylabel('normalized statistic','fontsize',16);
title('Sensitivity of Nonparametric Stats to Data Precision','fontsize',16);
grid on;

fname = 'figures/Figure6_wgewGaugeOnlySensitivity';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);



