clear all
close all
clc
addpath('tools');

% dimensions 
Ns = 50;
Nb = 9;
Nf = 10;

% number of bin sizes
B = linspace(5,100,Nb);
fracs = linspace(0,1,Nf+1); fracs(1) = [];

% read data
data = importdata('data/wgew_precip/gauge_data.txt');
date = data.data(:,1:3);
data = data.data(:,4:end);

% number of samples and gauges
[Nt,Ng] = size(data);

% truth
T = mean(data,2);

% random permutations
perms = zeros(Ns,3);

% loop through experiments
for s = 1:Ns

 % get a sample of gauges
 perms(s,:) = randperm(Ng,3);
 
 % create measurements 
 X = data(:,perms(s,1));
 Y = data(:,perms(s,2));
 Z = data(:,perms(s,3));

 % remove grandmas
 clear I
 I = find(any([X,Y,Z]'<=0));
 if ~isempty(I); 
  X(I) = []; 
  Y(I) = []; 
  Z(I) = []; 
 end;
 XX = X; YY = Y; ZZ = Z;
 N(s) = length(XX); 

 % loop through sample sizes
 for f = 1:Nf

  % select sample
  mm = length(XX);
  nn = round(fracs(f)*mm); 
  pp = randperm(mm,nn);
  X = XX(pp); Y = YY(pp); Z = ZZ(pp); 

  % loop through resolutions
  for b = 1:Nb

   % create bins at resolution
   Bt = linspace(min(T)-1e-6,max(T)+1e-6,B(b));
   Bx = linspace(min(X)-1e-6,max(X)+1e-6,B(b));
   By = linspace(min(Y)-1e-6,max(Y)+1e-6,B(b));
   Bz = linspace(min(Z)-1e-6,max(Z)+1e-6,B(b));

   % calculate nonlinear TC stats
   [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,Bx,By,Bz);

   % bound on total error
   NE(s,b,f,1,1) = 1-(Ixy+Ixz-Ixyz)/Hx;
   NE(s,b,f,2,1) = 1-(Ixy+Iyz-Ixyz)/Hy;
   NE(s,b,f,3,1) = 1-(Ixz+Iyz-Ixyz)/Hz;

   [s/Ns,f/Nf,b/Nb]

  end % bin resolution
 end % sample size
end % gauge triplets

% save results
save('results/wgew_convergece.mat');

%grab colors for ploting
figure(100); h = plot(randn(10,10));
colors = get(h,'Color');
close(100)

% legend names
for b = 1:Nb
 binNames(b) = strcat({'resolution: '},num2str(round(1/B(b)*1000)/1000),{' [m3/m3]'});
end

% plot convergence
figure(1); close(1); figure(1);
set(gcf,'color','w');%,'position',[1000 895 560 450]);

plotdata = squeeze(NE(:,:,:,:,1));
plotdat2 = cat(1,plotdata(:,:,:,1),plotdata(:,:,:,2));
plotdata = cat(1,plotdat2,plotdata(:,:,:,3));
mu = squeeze(mean(plotdata,1));

plot(fracs,mu','-o');%,'color',colors{b}); hold on; 
grid on;
if s == 5; legend(binNames','location','se'); end;
if s == 5; xlabel('number of data points','fontsize',16); end;
ylabel('H(X_i│T)/H(X_i)','fontsize',16); 


%% *******************************************
fname = 'figures/FigureX_WGEW_Convergence';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);


