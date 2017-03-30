clear all
close all
clc
addpath('tools');

% site names
siteNames = [{'Walnut Gulch'},{'Little Washita'},{'Fort Cobb'},{'Little River'},{'Reynolds Creek'}];
Nsites = length(siteNames);

% info bins
B = linspace(0,1,20);

% loop through sites
for s = 1:Nsites

 % screen report
 fprintf('Site %d of %d ... ',s,Nsites); tic;

 % ---- Loading Data --------------------------------------------
 % file name
 fname = strcat('data/smap_',num2str(s),'.txt');

 % load site data
 data = load(fname); 

 % replace missings with grandmas
 data(data<-9990) = 0/0;

 % find dates with missing SMAP, model, or averages
 data(any(isnan(data(:,4:7))'),:) = [];

 % find columns with more than half missing
 missing = zeros(size(data,2),1);
 for d = 8:size(data,2)
  I = find(isnan(data(:,d)));
  if length(I) > size(data,1)/2; missing(d) = 1; end;
 end 
 data(:,find(missing)) = [];

 % data dimensions
 [N(s),D(s)] = size(data);
 G(s) = D(s) - 7;

 % ---- Average Analysis ----------------------------------------
 % segregate data
 X = data(:,4); 
 Y = data(:,5); 
 Zs = data(:,6);
 Za = data(:,7);

 % linear TC
 [LEs(s,:),LIs(s,:)] = linear_triple_collocation(X,Y,Zs);
 [LEa(s,:),LIa(s,:)] = linear_triple_collocation(X,Y,Za);

 % nonlinear TC
 [NEs(s,:),NIs(s,:)] = nonlinear_triple_collocation(X,Y,Zs,B);
 [NEa(s,:),NIa(s,:)] = nonlinear_triple_collocation(X,Y,Za,B);

 % ---- Gauge Analysis ------------------------------------------
 % init storage
 LEg  = zeros(3,2,G(s)); 
 LIg  = zeros(3,2,G(s)); 
 NIg  = zeros(3,2,G(s)); 
 NEg  = zeros(3,2,G(s)); 
 IND  = zeros(3,G(s));

 remove = zeros(G(s),1);

 % loop through gauges
 for g = 1:G(s)

  % segregate data
  X = data(:,4); 
  Y = data(:,5); 
  Z = data(:,7+g); 
  T = data(:,6);
%  T = data(:,7);

  % remove missing data
  I = find(any(isnan([X,Y,Z]')));
  X(I) = []; Y(I) = []; Z(I) = []; T(I) = [];

  % linear TC 
  [LEg(:,1,g),LIg(:,1,g)] = linear_triple_collocation(X,Y,Z);

  % linear TC truths
  LEg(1,2,g) = cov(X-T)/cov(X);
  LEg(2,2,g) = cov(Y-T)/cov(Y);
  LEg(3,2,g) = cov(Z-T)/cov(Z);
  cc = corrcoef(X,T); LIg(1,2,g) = cc(2);
  cc = corrcoef(Y,T); LIg(2,2,g) = cc(2);
  cc = corrcoef(Z,T); LIg(3,2,g) = cc(2);

  % nonlinear tc routines
  [NEg(:,1,g),NIg(:,1,g)] = nonlinear_triple_collocation(X,Y,Z,B);

  % nonlinear tc truth
  [Ixt,Hx,Ht] = mutual_info(X,T,B,B);
  NIg(1,2,g) = Ixt/Hx;
  NEg(1,2,g) = 1-Ixt/Hx;

  [Ixt,Hx,Ht] = mutual_info(Y,T,B,B);
  NIg(2,2,g) = Ixt/Hx;
  NEg(2,2,g) = 1-Ixt/Hx;

  [Ixt,Hx,Ht] = mutual_info(Z,T,B,B);
  NIg(3,2,g) = Ixt/Hx;
  NEg(3,2,g) = 1-Ixt/Hx;

  % measure independence
  [Ixyt,Ixy,Ixt,Iyt,Hx,Hy,Ht] = mutual_info_3(X,Y,T,B,B,B);
  [Ixtz,Ixt,Ixz,Itz,Hx,Ht,Hz] = mutual_info_3(X,T,Z,B,B,B);
  [Ityz,Ity,Itz,Iyz,Ht,Hy,Hz] = mutual_info_3(T,Y,Z,B,B,B);
  IND(1,g) = (Ixy-Ixyt)/Hx;
  IND(2,g) = (Ixz-Ixtz)/Hy;
  IND(3,g) = (Iyz-Ityz)/Hz;

  % assert no errors
%  try
   gg = LEg(:,:,g); assert(~any(isnan(gg(:))));
   gg = LIg(:,:,g); assert(~any(isnan(gg(:))));
   gg = NEg(:,:,g); assert(~any(isnan(gg(:))));
   gg = NIg(:,:,g); assert(~any(isnan(gg(:))));
%  catch
%   remove(g) = 1;
%  end

 end

 % remove bad data
 LEg(:,:,find(remove)) = [];
 LIg(:,:,find(remove)) = [];
 NEg(:,:,find(remove)) = [];
 NIg(:,:,find(remove)) = [];
 IND(:,find(remove)) = [];

 % screen report
 t = toc; fprintf(' finished: time = %f \n',t);

 results(s).LEs = LEs; 
 results(s).LIs = LIs; 
 results(s).NEs = NEs; 
 results(s).NIs = NIs; 

 results(s).LEa = LEa; 
 results(s).LIa = LIa; 
 results(s).NEa = NEa; 
 results(s).NIa = NIa; 

 results(s).LEg = LEg; 
 results(s).LIg = LIg; 
 results(s).NEg = NEg; 
 results(s).NIg = NIg; 

 results(s).IND = IND;

end

% calc stats
ne_se_g = []; le_se_g = [];
ni_se_g = []; li_se_g = [];
ne_se_s = []; le_se_s = [];
ni_se_s = []; li_se_s = [];
ne_se_a = []; le_se_a = [];
ni_se_a = []; li_se_a = [];

for s = 1:Nsites

 ne_se_g = [ne_se_g,squeeze(results(s).NEg(:,1,:)-results(s).NEg(:,2,:)).^2];
 le_se_g = [le_se_g,squeeze(results(s).LEg(:,1,:)-results(s).LEg(:,2,:)).^2];
 ni_se_g = [ni_se_g,squeeze(results(s).NIg(:,1,:)-results(s).NIg(:,2,:)).^2];
 li_se_g = [li_se_g,squeeze(results(s).LIg(:,1,:)-results(s).LIg(:,2,:)).^2];

 ne_se_s = [ne_se_s,squeeze(results(s).NEg(:,1,:)-repmat(NEs(s,:)',[1,1,G(s)])).^2];
 le_se_s = [le_se_s,squeeze(results(s).LEg(:,1,:)-repmat(LEs(s,:)',[1,1,G(s)])).^2];
 ni_se_s = [ni_se_s,squeeze(results(s).NIg(:,1,:)-repmat(NIs(s,:)',[1,1,G(s)])).^2];
 li_se_s = [li_se_s,squeeze(results(s).LIg(:,1,:)-repmat(LIs(s,:)',[1,1,G(s)])).^2];

 ne_se_a = [ne_se_a,squeeze(results(s).NEg(:,1,:)-repmat(NEa(s,:)',[1,1,G(s)])).^2];
 le_se_a = [le_se_a,squeeze(results(s).LEg(:,1,:)-repmat(LEa(s,:)',[1,1,G(s)])).^2];
 ni_se_a = [ni_se_a,squeeze(results(s).NIg(:,1,:)-repmat(NIa(s,:)',[1,1,G(s)])).^2];
 li_se_a = [li_se_a,squeeze(results(s).LIg(:,1,:)-repmat(LIa(s,:)',[1,1,G(s)])).^2];

 [He(:,s),Pe(:,s)] = ttest2(log(squeeze(results(s).NEg(:,1,:)-results(s).NEg(:,2,:))'.^2),...
                            log(squeeze(results(s).LEg(:,1,:)-results(s).LEg(:,2,:))'.^2),'Alpha',0.1);
 De(:,s) = mean(squeeze(results(s).LEg(:,1,:)-results(s).LEg(:,2,:))'.^2)-...
           mean(squeeze(results(s).NEg(:,1,:)-results(s).NEg(:,2,:))'.^2);
 [Hi(:,s),Pi(:,s)] = ttest2(log(squeeze(results(s).NIg(:,1,:)-results(s).NIg(:,2,:))'.^2),...
                            log(squeeze(results(s).LIg(:,1,:)-results(s).LIg(:,2,:))'.^2),'Alpha',0.1);
 Di(:,s) = mean(squeeze(results(s).LIg(:,1,:)-results(s).LIg(:,2,:))'.^2)-...
           mean(squeeze(results(s).NIg(:,1,:)-results(s).NIg(:,2,:))'.^2);
end

ne_se_g(:,124) = []; le_se_g(:,124) = []; ni_se_g(:,124) = []; li_se_g(:,124) = [];
ne_se_g(:,15)  = []; le_se_g(:,15)  = []; ni_se_g(:,15)  = []; li_se_g(:,15)  = [];
ne_se_s(:,124) = []; le_se_s(:,124) = []; ni_se_s(:,124) = []; li_se_s(:,124) = [];
ne_se_s(:,15)  = []; le_se_s(:,15)  = []; ni_se_s(:,15)  = []; li_se_s(:,15)  = [];
ne_se_a(:,124) = []; le_se_a(:,124) = []; ni_se_a(:,124) = []; li_se_a(:,124) = [];
ne_se_a(:,15)  = []; le_se_a(:,15)  = []; ni_se_a(:,15)  = []; li_se_a(:,15)  = [];

[he,pe] = ttest2(log(ne_se_s'),log(le_se_s'),'Alpha',0.05);
[hi,pi] = ttest2(log(ni_se_s'),log(li_se_s'),'Alpha',0.05);
%[he,pe] = ttest2(ne_se',le_se','Alpha',0.10)
%[hi,pi] = ttest2(ni_se',li_se','Alpha',0.10)

[mean(le_se_g(1,:)),mean(ne_se_g(1,:)); ...
 mean(li_se_g(1,:)),mean(ni_se_g(1,:)); ...
 mean(le_se_g(2,:)),mean(ne_se_g(2,:)); ...
 mean(li_se_g(2,:)),mean(ni_se_g(2,:)); ...
 mean(le_se_g(3,:)),mean(ne_se_g(3,:)); ...   
 mean(li_se_g(3,:)),mean(ni_se_g(3,:))]   

[mean(le_se_a(1,:)),mean(ne_se_a(1,:)); ...
 mean(li_se_a(1,:)),mean(ni_se_a(1,:)); ...
 mean(le_se_a(2,:)),mean(ne_se_a(2,:)); ...
 mean(li_se_a(2,:)),mean(ni_se_a(2,:)); ...
 mean(le_se_a(3,:)),mean(ne_se_a(3,:)); ...   
 mean(li_se_a(3,:)),mean(ni_se_a(3,:))]   

[mean(le_se_s(1,:)),mean(ne_se_s(1,:)); ...
 mean(li_se_s(1,:)),mean(ni_se_s(1,:)); ...
 mean(le_se_s(2,:)),mean(ne_se_s(2,:)); ...
 mean(li_se_s(2,:)),mean(ni_se_s(2,:)); ...
 mean(le_se_s(3,:)),mean(ne_se_s(3,:)); ...   
 mean(li_se_s(3,:)),mean(ni_se_s(3,:))]   

[he;hi]'

%% *** Plot ***************************************************************
%grab colors for ploting
figure(100); h = plot(randn(10,10));
colors = get(h,'Color');
close(100);

% scatter plots
figure(1); close(1); figure(1);
set(gcf,'color','w','position',[36         309        1842         661]);

for s = 1:Nsites

 subplot(2,Nsites,s)
 plot(squeeze(results(s).LEg(1,2,:)),squeeze(results(s).LEg(1,1,:)),'o','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).LEg(2,2,:)),squeeze(results(s).LEg(2,1,:)),'o','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).LEg(3,2,:)),squeeze(results(s).LEg(3,1,:)),'o','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).LIg(1,2,:)),squeeze(results(s).LIg(1,1,:)),'+','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).LIg(2,2,:)),squeeze(results(s).LIg(2,1,:)),'+','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).LIg(3,2,:)),squeeze(results(s).LIg(3,1,:)),'+','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot([0,1],[0,1],'k--')
 grid on; 
 title(strcat(siteNames{s},': Linear'),'fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 axis([0,1,0,1]);
 set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);

 subplot(2,Nsites,Nsites+s)
 plot(squeeze(results(s).NEg(1,2,:)),squeeze(results(s).NEg(1,1,:)),'o','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).NEg(2,2,:)),squeeze(results(s).NEg(2,1,:)),'o','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).NEg(3,2,:)),squeeze(results(s).NEg(3,1,:)),'o','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).NIg(1,2,:)),squeeze(results(s).NIg(1,1,:)),'+','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).NIg(2,2,:)),squeeze(results(s).NIg(2,1,:)),'+','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(squeeze(results(s).NIg(3,2,:)),squeeze(results(s).NIg(3,1,:)),'+','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot([0,1],[0,1],'k--')
 grid on; 
 title(strcat(siteNames{s},': Nonlinear'),'fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 axis([0,1,0,1]);
 set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
 if s == 1; legend('error: SMAP','error: ECMWF','error: in situ','corr: SMAP','corr: ECMWF','corr: in situ','location','nw'); end

end

fname = strcat('figures/Figure8_SMAP_EMCWF_single');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

return
%% *** Plot ***************************************************************
% scatter plots
figure(2); close(2); figure(2);
set(gcf,'color','w','position',[36         309        1842         661]);

for s = 1:Nsites

 subplot(2,Nsites,s)
 plot(repmat(LEs(s,1),[G(s),1]),squeeze(results(s).LEg(1,1,:)),'o','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(repmat(LEs(s,2),[G(s),1]),squeeze(results(s).LEg(2,1,:)),'o','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(repmat(LEs(s,3),[G(s),1]),squeeze(results(s).LEg(3,1,:)),'o','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot(repmat(LIs(s,1),[G(s),1]),squeeze(results(s).LIg(1,1,:)),'+','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(repmat(LIs(s,2),[G(s),1]),squeeze(results(s).LIg(2,1,:)),'+','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(repmat(LIs(s,3),[G(s),1]),squeeze(results(s).LIg(3,1,:)),'+','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot([0,1],[0,1],'k--')
 grid on; 
 title(strcat(siteNames{s},': Linear'),'fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 axis([0,1,0,1]);
 set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);

 subplot(2,Nsites,Nsites+s)
 plot(repmat(NEs(s,1),[G(s),1]),squeeze(results(s).NEg(1,1,:)),'o','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(repmat(NEs(s,2),[G(s),1]),squeeze(results(s).NEg(2,1,:)),'o','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(repmat(NEs(s,3),[G(s),1]),squeeze(results(s).NEg(3,1,:)),'o','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot(repmat(NEa(s,1),[G(s),1]),squeeze(results(s).NIg(1,1,:)),'+','color',colors{1},'linewidth',2,'markersize',8); hold on;
 plot(repmat(NEa(s,2),[G(s),1]),squeeze(results(s).NIg(2,1,:)),'+','color',colors{2},'linewidth',2,'markersize',8); hold on;
 plot(repmat(NEa(s,3),[G(s),1]),squeeze(results(s).NIg(3,1,:)),'+','color',colors{3},'linewidth',2,'markersize',8); hold on;
 plot([0,1],[0,1],'k--')
 grid on; 
 title(strcat(siteNames{s},': Nonlinear'),'fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 axis([0,1,0,1]);
 set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
 if s == 1; legend('error: SMAP','error: ECMWF','error: in situ','corr: SMAP','corr: ECMWF','corr: in situ','location','nw'); end

end

fname = strcat('figures/Figure8_SMAP_EMCWF_core');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);



