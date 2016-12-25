clear all
close all
clc

% which sites to use
sites = [1601,1602,1603,1604,401];
Nsites = length(sites);

% loop through sites
for s = 1:Nsites

 % screen report
 fprintf('Site %d of %d ...',s,Nsites); tic;

 % --- distributed data -------------------------------
 % load distributed gauge data
 fname = strcat('smap_core/SoilMoistures_',num2str(sites(s)),'_60minutes.txt');
 gdata = load(fname);

 % change missings to grandmas 
 gdata(gdata<-9990) = 0./0;

 % remove data that we know is bad
 if s == 1; gdata(:,5+20) = []; end; % constant data from gauge
% if s == 5; gdata(:,5+17) = []; end; % constant data from gauge
% if s == 4; gdata(:,end-1) = []; end; % constant data from gauge

 % number of gauges
 [Ng,Dg] = size(gdata);
 Dg = Dg - 5; 

 % convert dates to day of year
 gdate = zeros(Ng,1);
 for g = 1:Ng
  gdate(g) = round(convertDate(gdata(g,1:5))*20)/20;
 end % convert dates

 % --- collocated data --------------------------------
 % load collocated data
 fname = strcat('core_data/',num2str(sites(s)),'3601.txt');
 cdata = load(fname);

 % change missings to grandmas 
 cdata(cdata<-9990) = 0./0;

 % number of timesteps
 Nc = size(cdata,1);

 % round day of year
 cdata(:,2) = round(cdata(:,2)*20)/20;

 % --- new data ---------------------------------------
 % init storage (year, doy, orbit, SMAP, ECMWF, WASM, mean, #Gauges)
 site = zeros(Nc,7+Dg)./0;

 % loop through times in averaged file
 for c = 1:Nc
  Iy = find(gdata(:,1) == cdata(c,1)); % year 
  Id = find(gdate(Iy)  == cdata(c,2)); % doy
  if ~isempty(Id)
   site(c,1:2) = cdata(c,1:2); % year, doy
   site(c,3)   = cdata(c,6);   % orbit mode
   site(c,4:5) = cdata(c,4:5); % SMAP, ECMWF
   site(c,6)   = cdata(c,3);   % areally-weighted gauge average
   site(c,7)   = nanmean(gdata(Iy(Id),6:end)); % arithmetic gauge average 
   site(c,8:end) = gdata(Iy(Id),6:end); % individual gauges
  end
 end % times in gdata file

 % save concatenated data
 fname = strcat('smap_',num2str(s),'.txt');
 save(fname,'site','-ascii');

 % screen report
 t = toc; fprintf('. finished: times = %d, gauges = %d, t = %f \n',Nc,Dg,t);

end % sites




