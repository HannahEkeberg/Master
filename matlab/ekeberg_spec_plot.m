clear all; clc;

% Load parsed matrices from ekeberg_spec_anal.m
load ni_mat.mat;
load fe_mat.mat;
load cu_mat.mat;
load ir_mat.mat;


% Matrix format:
%
% Each row corresponds to a key energy, from the ZZZ_key_energies matrices, 
%  or from your energy assignments list
%
% Column 1:  the key energy
% Column 2,10,18,..., 8.*(i-1) +2 : net gamma counts (efficiency and attenuation corrected)
% Column 3,11,19,..., 8.*(i-1) +3 : error (in %) in net gamma counts
% Column 4,12,20,..., 8.*(i-1) +4 : sample mass #
% Column 5,13,21,..., 8.*(i-1) +5 : spectrum count time (in Julian time)
% Column 6,14,22,..., 8.*(i-1) +6 : spectrum count length (in seconds of
% live time
% Column 7,15,23,..., 8.*(i-1) +7 : parent decay half-life (in s)
% Column 8,16,24,..., 8.*(i-1) +8 : gamma branching ratio (in %)
% Column 9,17,25,..., 8.*(i-1) +9 : uncertainty in gamma b.r. (in %)

ir_rhodrs = mean([55.9950   55.6010   55.6430   56.0000   55.1610   55.7310   56.6850   58.0300   56.6690   55.0650]);
cu_rhodrs = mean([22.3380   22.3250   22.3130   22.2840   22.4430   22.3960   22.3200   22.4010   22.4250   22.3140]);
ni_rhodrs = mean([22.7720   23.1180   22.3380   20.7040   21.7680   22.8610   23.0920   22.4090   21.7410   23.0930]);
fe_rhodrs = mean([20.030    20.017    19.948]);






% Set foils for analysis
data = ir_mat;
rhodrs = ir_rhodrs;
mu_attenuation = load('./ir_xcom.txt');
pp2 = pchip(mu_attenuation(:,1).*1e3,mu_attenuation(:,2));


num_spec = (size(data,2)-1)/8;



% % 33MeV Nickel foils:
rows_52mMn = [70];  %check, prob not a peak
rows_52Mn = [27,42,60,64,71]; %check
rows_54Mn = [35]; %check
rows_55Co = [3,13,15,18,30,34,41,46,57,63,67,69,84,98,105];  %check
rows_56Co = [29,45,55,59,66,82,90,91,92,99,101]; %check
rows_56Ni = [7,8,16,28,32,74]; %check
rows_57Co = [4,6,26]; %check
rows_57Ni = [26,12,25,40,49,58,68,80,81,86,87,95,104]; %check
rows_58Co = [31,38,77]; %check
rows_59Fe = [61]; %check
rows_60Cu = [14,17,23,36,43,44,48,52,62,83,85,88,89,93,97,100,102,103]; %check
rows_60mCo = [1]; %check
rows_61Cu = [2,9,11,19,20,22,24,33,37,39,47,50,51,54,56,72,76,79,94]; %check
rows_64Cu = [65];
rows_65Ni = [10,21,53,73,75,78];




% % 33MeV Copper foils:
rows_52Mn = [20,23,34,39]; %check
rows_56Co = [22,25,31,42,44,47]; %check
rows_57Co = [3,15];%check
rows_57Ni = [36];%check
rows_58Co = [21]; %check
rows_59Fe = [10, 26, 32]; %check
rows_60Co = [29,33]; %check
rows_61Co = [2];%check
rows_61Cu = [2,8,18,30]; %check
rows_62Zn = [1,4,6,7,9,11,12,14,16,17]; %check
rows_63Zn = [5,13,19,24,28,37,38,41,43,45,46]; %check
rows_64Cu = [35]; %check
rows_65Ni = [40]; %check
rows_65Zn = [27];%check


% % 33MeV Iron foils:
rows_48V = [24,26,37]; %check
rows_51Cr = [5]; %check
rows_51Mn = [12]; %check
rows_52mMn = [44]; %check
rows_52Mn = [6,11,19,23,35,39,45]; %check
rows_53Fe = [7]; %check
rows_54Mn = [17]; %check
rows_55Co = [1,8,14,22,38,42,43,55,16,50,54];
rows_56Co = [4,9,10,13,18,20,21,25,27,28,30,31,32,33,34,40,41,46,47,48,49,51,52,53];
rows_57Co = [2,3];
rows_58Co = [15];
rows_59Fe = [29,36];
%rows_55Co = [];


% % 33MeV Iridium foils:
rows_183Ta = [23,38]; %check
rows_186Re = [8,12]; %check
rows_186Ta = [20,58];
rows_187W = [54,68,73,9]; %check
rows_188Ir = [78,89,92,93,94];
rows_188mRe = [1,5];
rows_188Pt = [17,42,48,53];
rows_188Re = [77,82];
rows_189Ir = [3,21];
rows_189Pt = [2,6,22,31,35,60,75];
rows_189Re = [13,14,16,62];
rows_189W = [24,47,80]; 
rows_190Ir = [27,29,41,44,49,76,84]; 
rows_190mRe = [7,72];
rows_190Re = [50,79,88,90];
rows_191Pt = [4,15,25,37,39,45,51,59,70];
rows_192Ir = [11,18,19,26,30,33,34,40,46,52,56,57,63,66,67,81,85];
rows_193mPt = [10]; %the most important one!!!!
rows_194Ir = [28,32,64,69,71,83,86,87,91];
rows_194m2Ir = [36,43,55,61,65,74];




% Select rows to plot
% varToStr = @(x) inputname(1);
rows = rows_194m2Ir;
outName = './csv/194m2Ir';
% rows = 12;

% Gate on foil energy - energy = 0 searches on all energies
% energy = 0;
% energy = 913;

% loop over all energies
% for energy = 0:0   % Show all foils in one plot (not for analysis!)
%for energy = 228   % single foil
%for energy = 128:100:1028   % Just Nickel
%for energy = 129:100:1029   % Just Copper
%for energy = 126:100:326   % Just Iron
for energy = 177:100:1077   % Just Iridium

if energy==0
    masses = data(rows,8.*((1:num_spec)-1) +4);
    gammas = data(rows,8.*((1:num_spec)-1) +2);
    branching = data(rows,8.*((1:num_spec)-1) +8);
    err_gam = data(rows,8.*((1:num_spec)-1) +3);
    err_br  = data(rows,8.*((1:num_spec)-1) +9);
    
    ct = data(rows,8.*((1:num_spec)-1) +5);
    cl = data(rows,8.*((1:num_spec)-1) +6);
    lifetimes = data(rows,8.*((1:num_spec)-1) +7);
else
    [~,col] = find((data(rows,8.*((1:num_spec)-1) +4))==energy);
    % Extract unique columns
    col = unique(col');
    
%     8.*((col)-1) +4
    
    masses = data(rows,8.*(col-1) +4);
    gammas = data(rows,8.*(col-1) +2);
    branching = data(rows,8.*(col-1) +8);
    err_gam = data(rows,8.*(col-1) +3);
    err_br  = data(rows,8.*(col-1) +9);
    
    ct = data(rows,8.*(col-1) +5);
    cl = data(rows,8.*(col-1) +6);
    lifetimes = data(rows,8.*(col-1) +7);

end


lambdas = (log(2)./lifetimes);
I_gammas = branching.*.01;  %  in decimal, branching is in percent




% EoB Time
delta_ts = (ct - juliandate(datetime('26-Feb-2019 00:32:00')) ) .* 24;

% Note: 'gammas' contains the net peak area, corrected for efficiency and
% self-attenuation
activities = gammas.*lambdas./(1.* I_gammas.*(1-exp(-lambdas.*cl) ));

% errors = sqrt(  (err_gam.* gammas./ (30 .* branching)).^2 + ...
%          ( (gammas .* err_br) ./ (1 .* (branching.^2) )  ).^2);

err_counts = lambdas.*gammas.*(.01.*err_gam)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );
err_efficiency = lambdas.*gammas.*(.04)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );
err_I_gammas = lambdas.*gammas.*(err_br./branching)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );
err_live_time = lambdas.^2.*gammas.*(5)./( (I_gammas).*4 .*( sinh(lambdas.*cl.*0.5) ).^2  );
err_attenuation = lambdas.*gammas.*(0.05.*ppval(pp2,data(rows,1))).* (0.5e-3 .* rhodrs)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );  

errors = sqrt( err_counts.^2 + err_efficiency.^2 + err_I_gammas.^2 + err_live_time.^2 + err_attenuation.^2);
approximate_error = activities.*sqrt((0.01.*err_gam).^2 + (0.04).^2 + (err_br./branching).^2);

% Plot each time step as same color
% errorbar(delta_ts,activities,errors,'.')
% Plot each gamma line as same color
errorbar(delta_ts',activities',errors','.')
% Plot each time step as same color
% plot(delta_ts',activities','.')

% Dump to csv for python / gnuplot


outfile = [reshape(delta_ts,[1 numel(delta_ts)])' ,reshape(activities,[1 numel(activities)])', reshape(errors,[1 numel(errors)])'];
% Turn this line on to write files out!
csvwrite([outName '_' num2str(energy) '.dat'],outfile);

end




% write out liens of each plotted matrix into a single vector, cocatentate them into a ZZZ x 3 matreix, wrote to disk and fit in gnuplot for 
