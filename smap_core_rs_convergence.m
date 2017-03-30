clear all
close all
clc
addpath('tools');

% site names
siteNames = [{'Walnut Gulch'},{'Little Washita'},{'Fort Cobb'},{'Little River'},{'Reynolds Creek'}];
Nsites = length(siteNames);

% info bins
Nb = 3
B = [10,25,50];%round(linspace(5,50,Nb));

% data fractions
Nf = 10;
F = linspace(0.1,1,Nf);

% bootstrap
Nt = 10;

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

 % loop through sample sizes and bootstraps
 for f = 1:Nf
  for t = 1:Nt

   % segregate data
   P = randperm(N(s),round(F(f)*N(s)));
   X = data(P,4); 
   Y = data(P,5); 
   Zs = data(P,6);
   Za = data(P,7);

   % linear TC
   [LEs(s,f,:,t),LIs(s,f,:,t)] = linear_triple_collocation(X,Y,Zs);
   [LEa(s,f,:,t),LIa(s,f,:,t)] = linear_triple_collocation(X,Y,Za);

   % nonlinear TC
   for b = 1:Nb
    BB = linspace(0,1,B(b));
    [NEs(s,f,b,:,t),NIs(s,f,b,:,t)] = nonlinear_triple_collocation(X,Y,Zs,BB);
    [NEa(s,f,b,:,t),NIa(s,f,b,:,t)] = nonlinear_triple_collocation(X,Y,Za,BB);
   end

   % screen report
   t = toc; fprintf(' finished: time = %f \n',t);

  end
 end
end

%% *** Plot ***************************************************************
%grab colors for ploting
figure(100); h = plot(randn(10,10));
colors = get(h,'Color');
close(100);

NEa = mean(NEa,5);

for s = 1:Nsites
 der(s,:,1) = diff(NEa(s,:,2,1))./diff(F)./N(s) 
 der(s,:,2) = diff(NEa(s,:,2,2))./diff(F)./N(s) 
 der(s,:,3) = diff(NEa(s,:,2,3))./diff(F)./N(s) 
end

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[300,100,1500,1200]);
for s = 1:Nsites
 subplot(3,2,s)
% figure(s); close(s); figure(s);
% set(gcf,'color','w');%,'position',[1,9,2064,450]);
 h1 = plot(N(s)*F,squeeze(NEa(s,:,1,1)),'--s','linewidth',1,'color',colors{1}); hold on;
 h2 = plot(N(s)*F,squeeze(NEa(s,:,1,2)),'--s','linewidth',1,'color',colors{2}); hold on;
 h3 = plot(N(s)*F,squeeze(NEa(s,:,1,3)),'--s','linewidth',1,'color',colors{3}); hold on;
 h4 = plot(N(s)*F,squeeze(NEa(s,:,2,1)),'-o' ,'linewidth',3,'color',colors{1}); hold on;
 h5 = plot(N(s)*F,squeeze(NEa(s,:,2,2)),'-o' ,'linewidth',3,'color',colors{2}); hold on;
 h6 = plot(N(s)*F,squeeze(NEa(s,:,2,3)),'-o' ,'linewidth',3,'color',colors{3}); hold on;
 h7 = plot(N(s)*F,squeeze(NEa(s,:,3,1)),':^' ,'linewidth',1,'color',colors{1}); hold on;
 h8 = plot(N(s)*F,squeeze(NEa(s,:,3,2)),':^' ,'linewidth',1,'color',colors{2}); hold on;
 h9 = plot(N(s)*F,squeeze(NEa(s,:,3,3)),':^' ,'linewidth',1,'color',colors{3}); hold on;
 title(siteNames{s},'fontsize',16);
 grid on;
 if s == 5; legend('10% - SMAP','10% - ECMWF','10% - in situ',...
                   '4% - SMAP' ,'4% - ECMWF' ,'4% - in situ',...
                   '2% - SMAP' ,'2% - ECMWF' ,'2% - in situ',...
                   'location','se'); end;
% if s == 5; xlabel('number of data','fontsize',16); end;
 xlabel('sample size','fontsize',16); 
 ylabel('H(X_i│T)/H(X_i)','fontsize',16);
end

subplot(3,2,6)
plot(F(2:end),squeeze(der(:,:,1))','-o','linewidth',2,'color',colors{1}); hold on;
plot(F(2:end),squeeze(der(:,:,2))','-o','linewidth',2,'color',colors{2}); hold on;
plot(F(2:end),squeeze(der(:,:,3))','-o','linewidth',2,'color',colors{3}); hold on;
title('derivatives','fontsize',16);
xlabel('sample size fraction','fontsize',16)
ylabel('\Delta per sample','fontsize',16);
plot([F(2),F(end)],[0,0],'--k');


fname = strcat('figures/Figure8_convergence');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);



