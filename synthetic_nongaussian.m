clear all
close all
clc
restoredefaultpath
addpath(genpath(pwd))

%grab colors for ploting
figure(100); h = plot(randn(10,10));
colors = get(h,'Color');
close(100);

% scatter plots
figure(1); close(1); fig=figure(1);
set(gcf,'color','w','position',[3200,450,450,800]);

% dimensions of problem
Nt = 2e3;
Nk = 11;
Ns = 5;
Nb = 10;

Ntrials = 10;

% set varibale along each dimension
%alpha = linspace(0,2,Na);
%beta  = linspace(0,1,Nb);
K     = linspace(-2,0.5,Nk);
sigma = linspace(0.1,1,Ns);

% bins
B = linspace(6,100,Nb);

% set truth
T = randn(Nt,1);

%% **************************************************

% loop through trials
for t = 1:Ntrials

 % loop through experiments
 for k = 1:Nk
  for s = 1:Ns
 
   % generat errors
   pd1 = random('gev',K(k),sigma(s),0,Nt,1);
   pd2 = random('gev',K(k),sigma(s),0,Nt,1);
   pd3 = random('gev',K(k),sigma(s),0,Nt,1);
   %pd1 = makedist('Stable','alpha',alpha(a),'beta',beta(b),'gam',1,'delta',0);
   %pd2 = makedist('Stable','alpha',alpha(a),'beta',beta(b),'gam',1,'delta',0);
   %pd3 = makedist('Stable','alpha',alpha(a),'beta',beta(b),'gam',1,'delta',0);

   % create measurements 
   Xgev = T + pd1;
   Ygev = T + pd2; 
   Zgev = T + pd3;

   % generat errors
   pd1 = random('norm',0,std(pd1),Nt,1);
   pd2 = random('norm',0,std(pd2),Nt,1);
   pd3 = random('norm',0,std(pd3),Nt,1);

   % create measurements 
   Xgau = T + pd1;
   Ygau = T + pd2; 
   Zgau = T + pd3;

%   % calculate continuous linear TC stats
%   [LEgev(k,s,1,1,t),LEgev(k,s,2,1,t),LEgev(k,s,3,1,t),LIgev(k,s,1,1,t),LIgev(k,s,2,1,t),LIgev(k,s,3,1,t)] = ...
%     triple_collocation(Xgev,Ygev,Zgev);
%   [LEgau(k,s,1,1,t),LEgau(k,s,2,1,t),LEgau(k,s,3,1,t),LIgau(k,s,1,1,t),LIgau(k,s,2,1,t),LIgau(k,s,3,1,t)] = ...
%      triple_collocation(Xgau,Ygau,Zgau);
%
%   % calculate linear continuous truth
%   LEgev(k,s,1,2) = cov(Xgev-T)/cov(Xgev);
%   LEgev(k,s,2,2) = cov(Ygev-T)/cov(Ygev);
%   LEgev(k,s,3,2) = cov(Zgev-T)/cov(Zgev);
%   cc = corrcoef(Xgev,T); LIgev(k,s,1,2) = cc(2);
%   cc = corrcoef(Ygev,T); LIgev(k,s,2,2) = cc(2);
%   cc = corrcoef(Zgev,T); LIgev(k,s,3,2) = cc(2);
%
%   LEgau(k,s,1,2) = cov(Xgau-T)/cov(Xgau);
%   LEgau(k,s,2,2) = cov(Ygau-T)/cov(Ygau);
%   LEgau(k,s,3,2) = cov(Zgau-T)/cov(Zgau);
%   cc = corrcoef(Xgau,T); LIgau(k,s,1,2) = cc(2);
%   cc = corrcoef(Ygau,T); LIgau(k,s,2,2) = cc(2);
%   cc = corrcoef(Zgau,T); LIgau(k,s,3,2) = cc(2);

   % loop through resolutions
   for b = 1:Nb

    % create bins at resolution
    Bt = linspace(min(T)   -1e-6,max(T)   +1e-6,B(b));
    Bx = linspace(min(Xgev)-1e-6,max(Xgev)+1e-6,B(b));
    By = linspace(min(Ygev)-1e-6,max(Ygev)+1e-6,B(b));
    Bz = linspace(min(Zgev)-1e-6,max(Zgev)+1e-6,B(b));

    % calculate nonlinear TC stats
    [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(Xgev,Ygev,Zgev,Bx,By,Bz);

    % bound on total information
    NIgev(k,s,1,1,b,t) = (Ixy+Ixz-Ixyz)/Hx; 
    NIgev(k,s,2,1,b,t) = (Ixy+Iyz-Ixyz)/Hy; 
    NIgev(k,s,3,1,b,t) = (Ixz+Iyz-Ixyz)/Hz; 
    % bound on total error
    NEgev(k,s,1,1,b,t) = 1-(Ixy+Ixz-Ixyz)/Hx;
    NEgev(k,s,2,1,b,t) = 1-(Ixy+Iyz-Ixyz)/Hy;
    NEgev(k,s,3,1,b,t) = 1-(Ixz+Iyz-Ixyz)/Hz;
    % bound on missing information
    NMgev(k,s,1,1,b,t) = (Iyz-Ixyz)/Hx;
    NMgev(k,s,2,1,b,t) = (Ixz-Ixyz)/Hy;
    NMgev(k,s,3,1,b,t) = (Ixy-Ixyz)/Hz;

    % calculate true stats
    [Ixt,Hx,Ht] = mutual_info(Xgev,T,Bx,Bt);
    NIgev(k,s,1,2,b,t) = Ixt/Hx;
    NEgev(k,s,1,2,b,t) = 1-Ixt/Hx;
    NMgev(k,s,1,2,b,t) = (Ht-Ixt)/Hx;

    [Ixt,Hx,Ht] = mutual_info(Ygev,T,By,Bt);
    NIgev(k,s,2,2,b,t) = Ixt/Hx;
    NEgev(k,s,2,2,b,t) = 1-Ixt/Hx;
    NMgev(k,s,2,2,b,t) = (Ht-Ixt)/Hx;
  
    [Ixt,Hx,Ht] = mutual_info(Zgev,T,Bz,Bt);
    NIgev(k,s,3,2,b,t) = Ixt/Hx;
    NEgev(k,s,3,2,b,t) = 1-Ixt/Hx;
    NMgev(k,s,3,2,b,t) = (Ht-Ixt)/Hx;

   end % bin resolution

   % loop through resolutions
   for b = 1:Nb

    % create bins at resolution
    Bt = linspace(min(T)   -1e-6,max(T)   +1e-6,B(b));
    Bx = linspace(min(Xgau)-1e-6,max(Xgau)+1e-6,B(b));
    By = linspace(min(Ygau)-1e-6,max(Ygau)+1e-6,B(b));
    Bz = linspace(min(Zgau)-1e-6,max(Zgau)+1e-6,B(b));

    % calculate nonlinear TC stats
    [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(Xgau,Ygau,Zgau,Bx,By,Bz);

    % bound on total information
    NIgau(k,s,1,1,b,t) = (Ixy+Ixz-Ixyz)/Hx; 
    NIgau(k,s,2,1,b,t) = (Ixy+Iyz-Ixyz)/Hy; 
    NIgau(k,s,3,1,b,t) = (Ixz+Iyz-Ixyz)/Hz; 
    % bound on total error
    NEgau(k,s,1,1,b,t) = 1-(Ixy+Ixz-Ixyz)/Hx;
    NEgau(k,s,2,1,b,t) = 1-(Ixy+Iyz-Ixyz)/Hy;
    NEgau(k,s,3,1,b,t) = 1-(Ixz+Iyz-Ixyz)/Hz;
    % bound on missing information
    NMgau(k,s,1,1,b,t) = (Iyz-Ixyz)/Hx;
    NMgau(k,s,2,1,b,t) = (Ixz-Ixyz)/Hy;
    NMgau(k,s,3,1,b,t) = (Ixy-Ixyz)/Hz;

    % calculate true stats
    [Ixt,Hx,Ht] = mutual_info(Xgau,T,Bx,Bt);
    NIgau(k,s,1,2,b,t) = Ixt/Hx;
    NEgau(k,s,1,2,b,t) = 1-Ixt/Hx;
    NMgau(k,s,1,2,b,t) = (Ht-Ixt)/Hx;

    [Ixt,Hx,Ht] = mutual_info(Ygau,T,By,Bt);
    NIgau(k,s,2,2,b,t) = Ixt/Hx;
    NEgau(k,s,2,2,b,t) = 1-Ixt/Hx;
    NMgau(k,s,2,2,b,t) = (Ht-Ixt)/Hx;
  
    [Ixt,Hx,Ht] = mutual_info(Zgau,T,Bz,Bt);
    NIgau(k,s,3,2,b,t) = Ixt/Hx;
    NEgau(k,s,3,2,b,t) = 1-Ixt/Hx;
    NMgau(k,s,3,2,b,t) = (Ht-Ixt)/Hx;

   end % bin resolution

  [t/Ntrials,k/Nk,s/Ns]

  end
 end % error variance
end

NEgev = squeeze(mean(NEgev,6));
NEgau = squeeze(mean(NEgau,6));

%% **************************************************

for t = 1:Ntrials

 % loop through experiments
 for s = 1:Ns
 
  % generat errors
  pd1 = random('logn',0,sigma(s),Nt,1);
  pd2 = random('logn',0,sigma(s),Nt,1);
  pd3 = random('logn',0,sigma(s),Nt,1);

  % create measurements 
  Xlgn = T + pd1;
  Ylgn = T + pd2; 
  Zlgn = T + pd3;

  % generat errors
  pd1 = random('norm',0,std(pd1),Nt,1);
  pd2 = random('norm',0,std(pd2),Nt,1);
  pd3 = random('norm',0,std(pd3),Nt,1);

  % create measurements 
  Xga2 = T + pd1;
  Yga2 = T + pd2; 
  Zga2 = T + pd3;

%  % calculate continuous linear TC stats
%  [LElgn(s,1,1),LElgn(s,2,1),LElgn(s,3,1),LIlgn(s,1,1),LIlgn(s,2,1),LIlgn(s,3,1)] = ...
%     triple_collocation(Xlgn,Ylgn,Zlgn);
%  [LEga2(s,1,1),LEga2(s,2,1),LEga2(s,3,1),LIga2(s,1,1),LIga2(s,2,1),LIga2(s,3,1)] = ...
%     triple_collocation(Xga2,Yga2,Zga2);
%
%  % calculate linear continuous truth
%  LElgn(s,1,2) = cov(Xlgn-T)/cov(Xlgn);
%  LElgn(s,2,2) = cov(Ylgn-T)/cov(Ylgn);
%  LElgn(s,3,2) = cov(Zlgn-T)/cov(Zlgn);
%  cc = corrcoef(Xlgn,T); LIlgn(s,1,2) = cc(2);
%  cc = corrcoef(Ylgn,T); LIlgn(s,2,2) = cc(2);
%  cc = corrcoef(Zlgn,T); LIlgn(s,3,2) = cc(2);
%
%  LEga2(s,1,2) = cov(Xga2-T)/cov(Xga2);
%  LEga2(s,2,2) = cov(Yga2-T)/cov(Yga2);
%  LEga2(s,3,2) = cov(Zga2-T)/cov(Zga2);
%  cc = corrcoef(Xga2,T); LIga2(s,1,2) = cc(2);
%  cc = corrcoef(Yga2,T); LIga2(s,2,2) = cc(2);
%  cc = corrcoef(Zga2,T); LIga2(s,3,2) = cc(2);

  % loop through resolutions
  for b = 1:Nb

   % create bins at resolution
   Bt = linspace(min(T)   -1e-6,max(T)   +1e-6,B(b));
   Bx = linspace(min(Xlgn)-1e-6,max(Xlgn)+1e-6,B(b));
   By = linspace(min(Ylgn)-1e-6,max(Ylgn)+1e-6,B(b));
   Bz = linspace(min(Zlgn)-1e-6,max(Zlgn)+1e-6,B(b));

   % calculate nonlinear TC stats
   [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(Xlgn,Ylgn,Zlgn,Bx,By,Bz);

   % bound on total information
   NIlgn(s,1,1,b,t) = (Ixy+Ixz-Ixyz)/Hx; 
   NIlgn(s,2,1,b,t) = (Ixy+Iyz-Ixyz)/Hy; 
   NIlgn(s,3,1,b,t) = (Ixz+Iyz-Ixyz)/Hz; 
   % bound on total error
   NElgn(s,1,1,b,t) = 1-(Ixy+Ixz-Ixyz)/Hx;
   NElgn(s,2,1,b,t) = 1-(Ixy+Iyz-Ixyz)/Hy;
   NElgn(s,3,1,b,t) = 1-(Ixz+Iyz-Ixyz)/Hz;
   % bound on missing information
   NMlgn(s,1,1,b,t) = (Iyz-Ixyz)/Hx;
   NMlgn(s,2,1,b,t) = (Ixz-Ixyz)/Hy;
   NMlgn(s,3,1,b,t) = (Ixy-Ixyz)/Hz;

   % calculate true stats
   [Ixt,Hx,Ht] = mutual_info(Xlgn,T,Bx,Bt);
   NIlgn(s,1,2,b,t) = Ixt/Hx;
   NElgn(s,1,2,b,t) = 1-Ixt/Hx;
   NMlgn(s,1,2,b,t) = (Ht-Ixt)/Hx;

   [Ixt,Hx,Ht] = mutual_info(Ylgn,T,By,Bt);
   NIlgn(s,2,2,b,t) = Ixt/Hx;
   NElgn(s,2,2,b,t) = 1-Ixt/Hx;
   NMlgn(s,2,2,b,t) = (Ht-Ixt)/Hx;
 
   [Ixt,Hx,Ht] = mutual_info(Zlgn,T,Bz,Bt);
   NIlgn(s,3,2,b,t) = Ixt/Hx;
   NElgn(s,3,2,b,t) = 1-Ixt/Hx;
   NMlgn(s,3,2,b,t) = (Ht-Ixt)/Hx;

  end % bin resolution

  % loop through resolutions
  for b = 1:Nb

   % create bins at resolution
   Bt = linspace(min(T)   -1e-6,max(T)   +1e-6,B(b));
   Bx = linspace(min(Xga2)-1e-6,max(Xga2)+1e-6,B(b));
   By = linspace(min(Yga2)-1e-6,max(Yga2)+1e-6,B(b));
   Bz = linspace(min(Zga2)-1e-6,max(Zga2)+1e-6,B(b));

   % calculate nonlinear TC stats
   [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(Xga2,Yga2,Zga2,Bx,By,Bz);

   % bound on total information
   NIga2(s,1,1,b,t) = (Ixy+Ixz-Ixyz)/Hx; 
   NIga2(s,2,1,b,t) = (Ixy+Iyz-Ixyz)/Hy; 
   NIga2(s,3,1,b,t) = (Ixz+Iyz-Ixyz)/Hz; 
   % bound on total error
   NEga2(s,1,1,b,t) = 1-(Ixy+Ixz-Ixyz)/Hx;
   NEga2(s,2,1,b,t) = 1-(Ixy+Iyz-Ixyz)/Hy;
   NEga2(s,3,1,b,t) = 1-(Ixz+Iyz-Ixyz)/Hz;
   % bound on missing information
   NMga2(s,1,1,b,t) = (Iyz-Ixyz)/Hx;
   NMga2(s,2,1,b,t) = (Ixz-Ixyz)/Hy;
   NMga2(s,3,1,b,t) = (Ixy-Ixyz)/Hz;

   % calculate true stats
   [Ixt,Hx,Ht] = mutual_info(Xga2,T,Bx,Bt);
   NIga2(s,1,2,b,t) = Ixt/Hx;
   NEga2(s,1,2,b,t) = 1-Ixt/Hx;
   NMga2(s,1,2,b,t) = (Ht-Ixt)/Hx;

   [Ixt,Hx,Ht] = mutual_info(Yga2,T,By,Bt);
   NIga2(s,2,2,b,t) = Ixt/Hx;
   NEga2(s,2,2,b,t) = 1-Ixt/Hx;
   NMga2(s,2,2,b,t) = (Ht-Ixt)/Hx;
 
   [Ixt,Hx,Ht] = mutual_info(Zga2,T,Bz,Bt);
   NIga2(s,3,2,b,t) = Ixt/Hx;
   NEga2(s,3,2,b,t) = 1-Ixt/Hx;
   NMga2(s,3,2,b,t) = (Ht-Ixt)/Hx;

  end % bin resolution

 [t/Ntrials,k/Nk,s/Ns]

 end
end

NElgn = squeeze(mean(NElgn,6));
NEga2 = squeeze(mean(NEga2,6));

%plot results
subplot(2,1,1)

mu = squeeze(mean(abs(NEgev(:,:,1,2,:)-NEgev(:,:,1,1,:)),2));
plot(K,mu,'o-','linewidth',1,'color',colors{2}); hold on;

mu = squeeze(mean(abs(NEgev(:,:,1,2,:)-NEgau(:,:,1,2,:)),2));
plot(K,mu,'o-','linewidth',1,'color',colors{4}); hold on;

title('Generalized Extreme Value','fontsize',16);
xlabel('GEV shape parameter','fontsize',16);
ylabel('H(X_i│T)/H(X_i)','fontsize',16);
legend('estimator error','difference between GEV and Gaussian');

%plot results
subplot(2,1,1)

mu = squeeze(mean(abs(NElgn(:,:,1,2,:)-NElgn(:,:,1,1,:)),2));
plot(K,mu,'o-','linewidth',1,'color',colors{2}); hold on;

mu = squeeze(mean(abs(NElgn(:,:,1,2,:)-NEga2(:,:,1,2,:)),2));
plot(K,mu,'o-','linewidth',1,'color',colors{4}); hold on;

title('Log-Normal','fontsize',16);
xlabel('standard deviation','fontsize',16);
ylabel('H(X_i│T)/H(X_i)','fontsize',16);
legend('estimator error','difference between Log-Normal and Gaussian');

%% *******************************************
fname = 'figures/Figure5_LinearNonGaussian';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);


