clear all
close all
clc
addpath('tools');

%% *** Experimental Setup *************************************************

% site names
siteNames = [{'Reynolds Creek'},{'Walnut Gulch'},{'Little Washita'},{'Fort Cobb'},{'Little River'}];

% site IDs
sites = [401,1601:1604];
Nsites = length(sites);

% info bins
B = linspace(0,1,50);

% maximum number of triplets
Mtrips = 100;

%% *** Load Data **********************************************************

% loop through sites
for s = 1:Nsites

 % screen report
 fprintf('Loading data at site %d of %d ... ',s,Nsites);

 % file name
 fname = strcat('data/smap_core/SoilMoistures_',num2str(sites(s)),'_','60minutes.txt');

 % load site data
 data = load(fname); 

 % replace missings with grandmas
 data(data<0) = 0/0;

 % find dates with all missing
 data = data(any(~isnan(data(:,6:end))'),6:end);

 % find columns with more than half missing
 missing = zeros(size(data,2),1);
 for d = 1:size(data,2)
  I = find(isnan(data(:,d)));
  if length(I) > size(data,1)/2; missing(d) = 1; end;
 end 
 data(:,find(missing)) = [];

 % data dimensions
 [N(s),D(s)] = size(data);

 % store in cell
 store(s).data = data;

 % screen report
 fprintf('finished - #Samples = %d, #Data = %d \n',N(s),D(s));

end

%% *** Calculate Info Metrics *********************************************

% init storage
LE  = zeros(Nsites,Mtrips,3,2); 
LI  = zeros(Nsites,Mtrips,3,2); 
NI  = zeros(Nsites,Mtrips,3,2); 
NE  = zeros(Nsites,Mtrips,3,2); 
IND = zeros(Nsites,Mtrips,3); 

% loop through sites
for s = 1:Nsites

 % pull data
 data = store(s).data;

 % all possible combinations
 Ntrips(s) = nchoosek(D(s),3);
 for t = 1:Mtrips
  trips(t,:) = randperm(D(s),3);
  while any(trips(t,:)==20)
   trips(t,:) = randperm(D(s),3);
  end
 end

 % calculate truth
 truth = nanmean(data,2);

 % loop through trips
 for t = 1:Mtrips

  % segregate data
  X = data(:,trips(t,1)); 
  Y = data(:,trips(t,2)); 
  Z = data(:,trips(t,3)); 
  T = truth;

  % remove missing data
  I = find(any(isnan([X,Y,Z]')));
  X(I) = []; Y(I) = []; Z(I) = []; T(I) = [];

  % screen report
  fprintf('Site %d/%d - sample %d/%d/%d - Npoints = %d ... ',s,Nsites,t,Mtrips,Ntrips(s),length(X)); tic;

  % call linear tc routine
  [LE(s,t,1,1),LE(s,t,2,1),LE(s,t,3,1),LI(s,t,1,1),LI(s,t,2,1),LI(s,t,3,1)] = triple_collocation(X,Y,Z);

  % linear tc truths
  LE(s,t,1,2) = cov(X-T)/cov(X);
  LE(s,t,2,2) = cov(Y-T)/cov(Y);
  LE(s,t,3,2) = cov(Z-T)/cov(Z);
  cc = corrcoef(X,T); LI(s,t,1,2) = cc(2);
  cc = corrcoef(Y,T); LI(s,t,2,2) = cc(2);
  cc = corrcoef(Z,T); LI(s,t,3,2) = cc(2);

  % nonlinear tc routines
  [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,B,B,B);
  NI(s,t,1,1) = (Ixy+Ixz-Ixyz)/Hx;   
  NI(s,t,2,1) = (Ixy+Iyz-Ixyz)/Hy;   
  NI(s,t,3,1) = (Ixz+Iyz-Ixyz)/Hz;
  NE(s,t,1,1) = 1-(Ixy+Ixz-Ixyz)/Hx; 
  NE(s,t,2,1) = 1-(Ixy+Iyz-Ixyz)/Hy; 
  NE(s,t,3,1) = 1-(Ixz+Iyz-Ixyz)/Hz;

  % nonlinear tc truth
  [Ixt,Hx,Ht] = mutual_info(X,T,B,B);
  NI(s,t,1,2) = Ixt/Hx;
  NE(s,t,1,2) = 1-Ixt/Hx;

  [Ixt,Hx,Ht] = mutual_info(Y,T,B,B);
  NI(s,t,2,2) = Ixt/Hx;
  NE(s,t,2,2) = 1-Ixt/Hx;

  [Ixt,Hx,Ht] = mutual_info(Z,T,B,B);
  NI(s,t,3,2) = Ixt/Hx;
  NE(s,t,3,2) = 1-Ixt/Hx;

  % measure nonlinaer non-independence
  [Ixyt,Ixy,Ixt,Iyt,Hx,Hy,Ht] = mutual_info_3(X,Y,T,B,B,B);
  [Ixtz,Ixt,Ixt,Itz,Hx,Ht,Hz] = mutual_info_3(X,T,Z,B,B,B);
  [Ityz,Ity,Itz,Iyz,Ht,Hy,Hz] = mutual_info_3(T,Y,Z,B,B,B);
  IND(s,t,1) = (Ixy-Ixyt)/Ixy;
  IND(s,t,2) = (Ixz-Ixtz)/Ixz;
  IND(s,t,3) = (Iyz-Ityz)/Iyz;

  % measure linear non-independence


  % assert no errors
  assert(~any(isnan(LE(:))));
  assert(~any(isnan(LI(:))));
  assert(~any(isnan(NE(:))));
  assert(~any(isnan(NI(:))));

  % screen report
  t = toc; fprintf(' finished: time = %f \n',t);

 end
end

%% *** Save Results *******************************************************

save('results/smap_core_results.mat');

% calc stats
for s = 1:Nsites

 % loss functions (bias)
 net = NE(s,:,:,2); net = net(:);
 nee = NE(s,:,:,1); nee = nee(:);
 nit = NI(s,:,:,2); nit = nit(:);
 nie = NI(s,:,:,1); nie = nie(:);
 let = LE(s,:,:,2); let = let(:);
 lee = LE(s,:,:,1); lee = lee(:);
 lit = LI(s,:,:,2); lit = lit(:);
 lie = LI(s,:,:,1); lie = lie(:);

 ne_bias(s) = mean(abs(net-nee));
 ni_bias(s) = mean(abs(nit-nie));
 le_bias(s) = mean(abs(let-lee));
 li_bias(s) = mean(abs(lit-lie));

 % loss functions (variance)
 ne_var(s)  = mean((net-nee).^2);
 ni_var(s)  = mean((nit-nie).^2);
 le_var(s)  = mean((let-lee).^2);
 li_var(s)  = mean((lit-lie).^2);

 % loss functions (mse)
 ne_mse(s)  = ne_bias(s)^2+ne_var(s);
 ni_mse(s)  = ni_bias(s)^2+ni_var(s);
 le_mse(s)  = le_bias(s)^2+le_var(s);
 li_mse(s)  = li_bias(s)^2+li_var(s);

 % correlation
 cc = corrcoef(net,nee); ne_cor(s) = cc(2);
 cc = corrcoef(nit,nie); ni_cor(s) = cc(2);
 cc = corrcoef(let,lee); le_cor(s) = cc(2);
 cc = corrcoef(lit,lie); li_cor(s) = cc(2);

 % risk significance tests
 [he,pe(s)] = ttest(abs(net-nee),abs(let-lee),'Alpha',1e-3);
 [hi,pi(s)] = ttest(abs(nit-nie),abs(lit-lie),'Alpha',1e-3);
 hype(s) = he;
 hypi(s) = hi;

 [mean(abs(net-nee)-abs(let-lee)),mean(abs(nit-nie)-abs(lit-lie))]

end

% screen repot
[hype;pe]
[hypi;pi]

% screen report
error_bias = [ne_bias',le_bias']
info_bias  = [ni_bias',li_bias']
error_var  = [ne_var' ,le_var' ]
info_var   = [ni_var' ,li_var' ]
error_mse  = [ne_mse' ,le_mse' ]
info_mse   = [ni_mse' ,li_mse' ]
error_cor  = [ne_cor' ,le_cor' ]
info_cor   = [ni_cor' ,li_cor' ]

% concatenate
for s = 1:Nsites
 stats(1,:,s) = [le_bias,ne_bias]; 
 stats(2,:,s) = [le_var,ne_var]; 
 stats(3,:,s) = [le_mse,ne_mse]; 
end
stats

% signal to noise ratios
s1 =  LI(:,:,:,2);              s1 = mean(s1(:));
s2 = (LI(:,:,:,2)-LI(:,:,:,1)); s2 = std(s2(:));
s2/s1

s1 = LE(:,:,:,2);               s1 = mean(s1(:));
s2 = (LE(:,:,:,2)-LE(:,:,:,1)); s2 = std(s2(:));
s2/s1

s1 = NI(:,:,:,2);               s1 = mean(s1(:));
s2 = (NI(:,:,:,2)-NI(:,:,:,1)); s2 = std(s2(:));
s2/s1

s1 = NE(:,:,:,2);               s1 = mean(s1(:));
s2 = (NE(:,:,:,2)-NE(:,:,:,1)); s2 = std(s2(:));
s2/s1

%% *** Plot ***************************************************************

%grab colors for ploting
figure(100); h = plot(randn(10,10));
colors = get(h,'Color');

% scatter plots
figure(1); close(1); figure(1);
set(gcf,'color','w','position',[36         309        1842         661]);

for s = 1:Nsites

 subplot(2,5,s)
 h1 = plot(squeeze(LE(s,:,:,2)),squeeze(LE(s,:,:,1)),'o','color',colors{4}); hold on;
 h2 = plot(squeeze(LI(s,:,:,2)),squeeze(LI(s,:,:,1)),'+','color',colors{2}); hold on;
 plot([0,1],[0,1],'k--')
 grid on; 
 title(strcat(siteNames{s},': Linear'),'fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 if s == 1; legend([h1(1),h2(1)],'total error','correlation','location','nw'); end
 axis([0,1,0,1]);

 subplot(2,5,5+s)
 h1 = plot(squeeze(NE(s,:,:,2)),squeeze(NE(s,:,:,1)),'o','color',colors{4}); hold on;
 h2 = plot(squeeze(NI(s,:,:,2)),squeeze(NI(s,:,:,1)),'+','color',colors{2}); hold on;
 plot([0,1],[0,1],'k--')
 grid on; 
 title(strcat(siteNames{s},': Nonlinear'),'fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 axis([0,1,0,1]);

end

fname = strcat('figures/Figure6_',siteNames{s},'_SMAPCoreGauges');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);


% scatter plots
figure(2); close(2); figure(2);
set(gcf,'color','w');

for s = 1:Nsites

 subplot(2,1,1)
 h1 = plot(squeeze(LE(s,:,:,2)),squeeze(LE(s,:,:,1)),'+','color',colors{s}); hold on;
 h2 = plot(squeeze(LI(s,:,:,2)),squeeze(LI(s,:,:,1)),'o','color',colors{s}); hold on;
 plot([0,1],[0,1],'k--')
 grid on;
 title('Linear TC','fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 if s == 1; legend([h1(1),h2(1)],'total error','correlation','location','nw'); end
 axis([0,1,0,1]);

 subplot(2,1,2)
 h1 = plot(squeeze(NE(s,:,:,2)),squeeze(NE(s,:,:,1)),'+','color',colors{s}); hold on;
 h2 = plot(squeeze(NI(s,:,:,2)),squeeze(NI(s,:,:,1)),'o','color',colors{s}); hold on;
 plot([0,1],[0,1],'k--')
 grid on;
 title('Nonlinear TC','fontsize',16);
 xlabel('true statistic','fontsize',18);
 ylabel('estimated statistic','fontsize',18);
 axis([0,1,0,1]);

 fname = strcat('figures/Figure6b_',siteNames{s},'_SMAPCoreGauges');
 img = getframe(gcf);
 imwrite(img.cdata, [fname, '.png']);

end



