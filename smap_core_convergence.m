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

% maximum number of triplets
Mtrips = 100;

% number of bin sizes
Nb = 10;
bins = linspace(10,100,Nb);

% number of sample sizes
Nn = 10;
fracs = linspace(0,1,Nn+1); fracs(1) = [];

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
NI  = zeros(Nsites,Mtrips,Nb,Nn,3,2); 
NE  = zeros(Nsites,Mtrips,Nb,Nn,3,2); 
IND = zeros(Nsites,Mtrips,Nb,Nn,3); 

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

  % loop through sample sizes
  for n = 1:Nn

   X = data(:,trips(t,1)); 
   Y = data(:,trips(t,2)); 
   Z = data(:,trips(t,3)); 
   T = truth;

   % randomize sample
   mm = length(X);
   nn = round(fracs(n)*mm); 
   pp = randperm(mm,nn);
   X = X(pp); Y = Y(pp); Z = Z(pp); T = T(pp);

   % segregate data
   % remove missing data
   I = find(any(isnan([X,Y,Z]')));
   X(I) = []; Y(I) = []; Z(I) = []; T(I) = [];

   % call linear tc routine
   [LE(s,t,1,1),LE(s,t,2,1),LE(s,t,3,1),LI(s,t,1,1),LI(s,t,2,1),LI(s,t,3,1)] = triple_collocation(X,Y,Z);

   % linear tc truths
   LE(s,t,1,2) = cov(X-T)/cov(X);
   LE(s,t,2,2) = cov(Y-T)/cov(Y);
   LE(s,t,3,2) = cov(Z-T)/cov(Z);
   cc = corrcoef(X,T); LI(s,t,1,2) = cc(2);
   cc = corrcoef(Y,T); LI(s,t,2,2) = cc(2);
   cc = corrcoef(Z,T); LI(s,t,3,2) = cc(2);

   % loop through bin sizes
   for b = 1:Nb
   
    % screen report
    fprintf('Site %d/%d - sample %d/%d/%d - Npoints = %d/%d - Bin %d/%d - Fraction %d/%d ... ',s,Nsites,t,Mtrips,Ntrips(s),length(X),mm,b,Nb,n,Nn); tic;

    % info bins
    B = linspace(0,1,bins(b));

    % nonlinear tc routines
    [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,B,B,B);
    NI(s,t,b,n,1,1) = (Ixy+Ixz-Ixyz)/Hx;   
    NI(s,t,b,n,2,1) = (Ixy+Iyz-Ixyz)/Hy;   
    NI(s,t,b,n,3,1) = (Ixz+Iyz-Ixyz)/Hz;
    NE(s,t,b,n,1,1) = 1-(Ixy+Ixz-Ixyz)/Hx; 
    NE(s,t,b,n,2,1) = 1-(Ixy+Iyz-Ixyz)/Hy; 
    NE(s,t,b,n,3,1) = 1-(Ixz+Iyz-Ixyz)/Hz;

    % nonlinear tc truth
    [Ixt,Hx,Ht] = mutual_info(X,T,B,B);
    NI(s,t,b,n,1,2) = Ixt/Hx;
    NE(s,t,b,n,1,2) = 1-Ixt/Hx;

    [Ixt,Hx,Ht] = mutual_info(Y,T,B,B);
    NI(s,t,b,n,2,2) = Ixt/Hx;
    NE(s,t,b,n,2,2) = 1-Ixt/Hx;

    [Ixt,Hx,Ht] = mutual_info(Z,T,B,B);
    NI(s,t,b,n,3,2) = Ixt/Hx;
    NE(s,t,b,n,3,2) = 1-Ixt/Hx;

    % measure nonlinaer non-independence
    [Ixyt,Ixy,Ixt,Iyt,Hx,Hy,Ht] = mutual_info_3(X,Y,T,B,B,B);
    [Ixtz,Ixt,Ixt,Itz,Hx,Ht,Hz] = mutual_info_3(X,T,Z,B,B,B);
    [Ityz,Ity,Itz,Iyz,Ht,Hy,Hz] = mutual_info_3(T,Y,Z,B,B,B);
    IND(s,t,b,n,1) = (Ixy-Ixyt)/Ixy;
    IND(s,t,b,n,2) = (Ixz-Ixtz)/Ixz;
    IND(s,t,b,n,3) = (Iyz-Ityz)/Iyz;

    % measure linear non-independence


    % assert no errors
    assert(~any(isnan(LE(:))));
    assert(~any(isnan(LI(:))));
    assert(~any(isnan(NE(:))));
    assert(~any(isnan(NI(:))));

    % screen report
    time = toc; fprintf(' finished: time = %f \n',time);

   end % bins
  end % trips
 end % sample size 
end % site

%% *** Save Results *******************************************************

save('results/smap_core_convergence_results.mat');

%% *** Plot ***************************************************************

%grab colors for ploting
figure(100); h = plot(randn(10,10));
colors = get(h,'Color');
close(100)

%[XX,YY] = meshgrid(bins,fracs);

for b = 1:Nb
 binNames(b) = strcat({'resolution: '},num2str(round(1/bins(b)*1000)/1000),{' [m3/m3]'});
end

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[1000 895 560 450]);

for s = 1:Nsites

 %NE(s,t,b,n,1,1) = 1-(Ixy+Ixz-Ixyz)/Hx; 
 plotdata = squeeze(NE(s,:,:,:,:,1));
 plotdat2 = cat(1,plotdata(:,:,:,1),plotdata(:,:,:,2));
 plotdata = cat(1,plotdat2,plotdata(:,:,:,3));
 mu = squeeze(mean(plotdata,1));
 sg = squeeze(std(plotdata,[],1));

 subplot(Nsites,1,s);
 plot(N(s)*fracs,mu','-o');%,'color',colors{b}); hold on; 
 grid on;
 title(siteNames{s},'fontsize',16);
 if s == 5; legend(binNames','location','se'); end;
 if s == 5; xlabel('number of data points','fontsize',16); end;
 ylabel('H(X_i│T)/H(X_i)','fontsize',16); 

end


%% *******************************************
fname = 'figures/FigureX_SMAPCore_Convergence';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);









