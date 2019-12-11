clear all; clc;
fclose('all');

load efficiency_matrix.mat

load fe_key_energies.mat
load cu_key_energies.mat
load ni_key_energies.mat
load ir_key_energies.mat
% 
% % Gamma line half-lives and intensities
load fe_glines.mat
load cu_glines.mat
load ni_glines.mat
load ir_glines.mat



% Load file names


%
% Change element prefix for energy identification, and pre-allocate size
%
%
fe_fn = textscan(fopen('fe_fnames.txt'),'%s');   
cu_fn = textscan(fopen('cu_fnames.txt'),'%s');   
ir_fn = textscan(fopen('ir_fnames.txt'),'%s');   
ni_fn = textscan(fopen('ni_fnames.txt'),'%s');   


% Pre-allocate size of arrays
% 33 MeV
fe_energies = zeros(2000,1);
cu_energies = zeros(2000,1);
ni_energies = zeros(2000,1);
ir_energies = zeros(2000,1);

%ir_glines = zeros(2000,3);
%ir_key_energies = zeros(2000,1);
%ir_mat = zeros(2000,1225);


ir_rhodrs = [55.9950   55.6010   55.6430   56.0000   55.1610   55.7310   56.6850   58.0300   56.6690   55.0650];  %mass densities
cu_rhodrs = [22.3380   22.3250   22.3130   22.2840   22.4430   22.3960   22.3200   22.4010   22.4250   22.3140];
ni_rhodrs = [22.7720   23.1180   22.3380   20.7040   21.7680   22.8610   23.0920   22.4090   21.7410   23.0930];
fe_rhodrs = [20.030    20.017    19.948];



% foil selection !
fn = ir_fn;   %fn=filename
energies = ir_energies;
rhodrs = ir_rhodrs;
mu_attenuation = load('../ir_xcom.txt');
% foil selection !


%% 
% Fit attenuation data to smooth function
%
% pp = spline(mu_attenuation(:,1).*1e3,mu_attenuation(:,2));
pp2 = pchip(mu_attenuation(:,1).*1e3,mu_attenuation(:,2));
xp=logspace(0,4,10000);
% yp=ppval(pp,xp);
yp2=ppval(pp2,xp);
loglog(mu_attenuation(:,1).*1e3,mu_attenuation(:,2),'--',xp,yp2);
xlabel('Gamma-Ray Energy (keV)')
ylabel('\mu (cm^2/g)')
clf;



%
% End of selection
%



%% Plotting extracted unique peak energies



% Changing indices for vector insertion
a = 1;  %global indices for final array we are making
b = 1;

for i=1:length(fn{1})
% for i=1:1
    % for i=1:3
    fname = fn{1}{i,1};
    fid = fopen(fname);
    
%     C2 = textscan(fid,'%d %f %f %f %f %d %f %f %f %f %f %s','headerlines',30);
    C2 = textscan(fid,'%f %f %f %f %f %d %f %f %f %d','headerlines',37);
    b = a + length(C2{1,1}) - 1;
    
    energies(a:b) =  C2{1,1};
    a = a + length(C2{1,1});
    
    fclose(fid);
end

energies = sort(energies);

%
% Be careful - veto  peaks
%
% energies(913) = [];
% energies(1438) = [];
% energies(773) = [];
%
%  End being careful
%

xp = 1:length(energies);
% plot(1:length(nb_energies),nb_energies,'.')
%
% Set tolerance (in keV) for peak uniqueness
%
tol = 0.75;

[C,IA] = uniquetol(energies, tol , 'DataScale', 1,'OutputAllIndices',true);


hold on
for k = 1:length(IA)
    plot(xp(IA{k}), energies(IA{k}), '.')
    meanAi = mean(energies(IA{k},:));
    C(k) = meanAi;
    plot(find((energies-meanAi)>=0,1), meanAi, 'xr')
end

xlabel('Peak Index')
ylabel('Gamma-Ray Energy (keV)')
title('Ir')
unique_energies = C;

% close all;


% nb_fn = string(nb_fn{1})
% Extract name i and remove .dat extension
% nb_fn{1}{i,1}(1:length(nb_fn{1}{i,1})-4)







%%
% % 
% not using this block
% % unique shelves:  1     3     5    10    14    15    18    22    30    40    60
% % 
% 
% eff_30 = @(x)(exp(effi_mat(1,2).*(log(x)).^2 + effi_mat(1,3).*log(x) + effi_mat(1,4)))  ;
% eff_1 = @(x)(exp(effi_mat(2,2).*(log(x)).^2 + effi_mat(2,3).*log(x) + effi_mat(2,4)))  ;
% eff_5 = @(x)(exp(effi_mat(3,2).*(log(x)).^2 + effi_mat(3,3).*log(x) + effi_mat(3,4)))  ;
% eff_40 = @(x)(exp(effi_mat(4,2).*(log(x)).^2 + effi_mat(4,3).*log(x) + effi_mat(4,4)))  ;
% eff_10 = @(x)(exp(effi_mat(5,2).*(log(x)).^2 + effi_mat(5,3).*log(x) + effi_mat(5,4)))  ;
% eff_18 = @(x)(exp(effi_mat(6,2).*(log(x)).^2 + effi_mat(6,3).*log(x) + effi_mat(6,4)))  ;
% eff_3 = @(x)(exp(effi_mat(7,2).*(log(x)).^2 + effi_mat(7,3).*log(x) + effi_mat(7,4)))  ;
% eff_10n = @(x)(exp(effi_mat(8,2).*(log(x)).^2 + effi_mat(8,3).*log(x) + effi_mat(8,4)))  ;
% eff_12 = @(x)(exp(effi_mat(9,2).*(log(x)).^2 + effi_mat(9,3).*log(x) + effi_mat(9,4)))  ;
% eff_14 = @(x)(exp(effi_mat(10,2).*(log(x)).^2 + effi_mat(10,3).*log(x) + effi_mat(10,4)))  ;
% eff_15 = @(x)(exp(effi_mat(11,2).*(log(x)).^2 + effi_mat(11,3).*log(x) + effi_mat(11,4)))  ;
% eff_22 = @(x)(exp(effi_mat(12,2).*(log(x)).^2 + effi_mat(12,3).*log(x) + effi_mat(12,4)))  ;
% eff_60 = @(x)(exp(effi_mat(13,2).*(log(x)).^2 + effi_mat(13,3).*log(x) + effi_mat(13,4)))  ;
% 
% beta=[ 0.1613   -2.8769    7.6767];
% eff_3n = @(x)(exp(beta(1).*(log(x)).^2 + beta(2).*log(x) + beta(3)))  ;
% 
% xp = 62:1:2500;
% 
% % semilogy(xp, eff_1(xp),xp, eff_5(xp), xp, sqrt(eff_1(xp).*eff_5(xp)), xp, eff_3n(xp))
% % semilogy( xp, sqrt(eff_1(xp).*eff_5(xp)), xp, eff_3n(xp))
% 
% % plot(xp, sqrt(eff_1(xp).*eff_5(xp)) -  eff_3n(xp))
% % plot(xp, (sqrt(eff_1(xp).*eff_5(xp)) -  eff_3n(xp))./eff_3n(xp))
% 
% semilogy(xp, eff_1(xp),xp, eff_5(xp),xp, eff_10(xp), xp, eff_10n(xp),xp, eff_14(xp),xp, eff_15(xp),xp, eff_18(xp),...
%     xp, eff_22(xp),xp, eff_30(xp),xp, eff_40(xp),xp, eff_60(xp), xp, sqrt(eff_1(xp).*eff_5(xp)), xp, eff_3n(xp))
% 
% ylim([6E-5 3E-1])

%%

clc;

% Convert to array of strings
fe_fn = string(fe_fn{1});
cu_fn = string(cu_fn{1});
ni_fn = string(ni_fn{1});
ir_fn = string(ir_fn{1});

% Select foils to analyze
fn = ir_fn;
key_energies = ir_key_energies;
glines = ir_glines;



data_mat = zeros(length(key_energies), 1+ (8 * length(fn)));
data_mat(:,1) = key_energies;

shelves = zeros(1,length(fn));
detectors = shelves;


for i=1:length(fn)
    
    fname = char(fn(i));
    fid = fopen(fname);
    
    % Get raw text for header regex extraction
    raw_str = fileread(char(fn(i)));
    
    % Parse column data to cell structure
    %   '37' refers to the number of header lines in the report file 
    C3 = textscan(fid,'%f %f %f %f %f %d %f %f %f %d','headerlines',37);
    fclose(fid);
    
    % Get shelf position, for efficiency correction
    shelf = regexp(raw_str, 'Shelf:\s+(\d+)[on]?', 'tokens');
    shelves(i) = cell2mat(cellfun(@(x) str2double(x{:}), shelf, 'UniformOutput', false));
    shelf = shelves(i);
    
    % Get detector ID, for efficiency correction
    detector = regexp(raw_str, 'Detector:\s+(\d+)[on]?', 'tokens');
    detectors(i) = cell2mat(cellfun(@(x) str2double(x{:}), detector, 'UniformOutput', false));
    detector = detectors(i);
    
%     Select efficiency curve based on shelf
    if detector==1
        if shelf==10
            effcal = @(x)(efficiency_matrix(1,1).*exp(-efficiency_matrix(1,2).*x.^efficiency_matrix(1,3)) .* (1-exp(-efficiency_matrix(1,4).*x.^efficiency_matrix(1,5)))  );
        elseif shelf==30
            effcal = @(x)(efficiency_matrix(2,1).*exp(-efficiency_matrix(2,2).*x.^efficiency_matrix(2,3)) .* (1-exp(-efficiency_matrix(2,4).*x.^efficiency_matrix(2,5)))  );
        elseif shelf==50
            effcal = @(x)((32/51.3).^2.*efficiency_matrix(2,1).*exp(-efficiency_matrix(2,2).*x.^efficiency_matrix(2,3)) .* (1-exp(-efficiency_matrix(2,4).*x.^efficiency_matrix(2,5)))  );
        end
    elseif detector==2
        if shelf==10
            effcal = @(x)(efficiency_matrix(3,1).*exp(-efficiency_matrix(3,2).*x.^efficiency_matrix(3,3)) .* (1-exp(-efficiency_matrix(3,4).*x.^efficiency_matrix(3,5)))  );
        elseif shelf==30
            effcal = @(x)(efficiency_matrix(4,1).*exp(-efficiency_matrix(4,2).*x.^efficiency_matrix(4,3)) .* (1-exp(-efficiency_matrix(4,4).*x.^efficiency_matrix(4,5)))  );
        elseif shelf==45
            effcal = @(x)((34.8/45.4).^2.*efficiency_matrix(4,1).*exp(-efficiency_matrix(4,2).*x.^efficiency_matrix(4,3)) .* (1-exp(-efficiency_matrix(4,4).*x.^efficiency_matrix(4,5)))  );
        end
    elseif detector==3
        if shelf==53
            effcal = @(x)(efficiency_matrix(5,1).*exp(-efficiency_matrix(5,2).*x.^efficiency_matrix(5,3)) .* (1-exp(-efficiency_matrix(5,4).*x.^efficiency_matrix(5,5)))  );
        end
    elseif detector==4
        if shelf==32
            effcal = @(x)(efficiency_matrix(6,1).*exp(-efficiency_matrix(6,2).*x.^efficiency_matrix(6,3)) .* (1-exp(-efficiency_matrix(6,4).*x.^efficiency_matrix(6,5)))  );
        end
    elseif detector==5
        if shelf==40
            effcal = @(x)(efficiency_matrix(7,1).*exp(-efficiency_matrix(7,2).*x.^efficiency_matrix(7,3)) .* (1-exp(-efficiency_matrix(7,4).*x.^efficiency_matrix(7,5)))  );
        end
    elseif detector==6
        if shelf==25
            effcal = @(x)(efficiency_matrix(8,1).*exp(-efficiency_matrix(8,2).*x.^efficiency_matrix(8,3)) .* (1-exp(-efficiency_matrix(8,4).*x.^efficiency_matrix(8,5)))  );
        end
    elseif detector==7
        if shelf==5
            effcal = @(x)(efficiency_matrix(9,1).*exp(-efficiency_matrix(9,2).*x.^efficiency_matrix(9,3)) .* (1-exp(-efficiency_matrix(9,4).*x.^efficiency_matrix(9,5)))  );
        elseif shelf==10
            effcal = @(x)(efficiency_matrix(10,1).*exp(-efficiency_matrix(10,2).*x.^efficiency_matrix(10,3)) .* (1-exp(-efficiency_matrix(10,4).*x.^efficiency_matrix(10,5)))  );
        elseif shelf==15
            effcal = @(x)(efficiency_matrix(11,1).*exp(-efficiency_matrix(11,2).*x.^efficiency_matrix(11,3)) .* (1-exp(-efficiency_matrix(11,4).*x.^efficiency_matrix(11,5)))  );
        elseif shelf==18
            effcal = @(x)(efficiency_matrix(12,1).*exp(-efficiency_matrix(12,2).*x.^efficiency_matrix(12,3)) .* (1-exp(-efficiency_matrix(12,4).*x.^efficiency_matrix(12,5)))  );    
        elseif shelf==22
            effcal = @(x)(efficiency_matrix(13,1).*exp(-efficiency_matrix(13,2).*x.^efficiency_matrix(13,3)) .* (1-exp(-efficiency_matrix(13,4).*x.^efficiency_matrix(13,5)))  );
        elseif shelf==30
            effcal = @(x)(efficiency_matrix(14,1).*exp(-efficiency_matrix(14,2).*x.^efficiency_matrix(14,3)) .* (1-exp(-efficiency_matrix(14,4).*x.^efficiency_matrix(14,5)))  );
        elseif shelf==40
            effcal = @(x)(efficiency_matrix(15,1).*exp(-efficiency_matrix(15,2).*x.^efficiency_matrix(15,3)) .* (1-exp(-efficiency_matrix(15,4).*x.^efficiency_matrix(15,5)))  );
        elseif shelf==50
            effcal = @(x)(efficiency_matrix(16,1).*exp(-efficiency_matrix(16,2).*x.^efficiency_matrix(16,3)) .* (1-exp(-efficiency_matrix(16,4).*x.^efficiency_matrix(16,5)))  );
        elseif shelf==60
            effcal = @(x)(efficiency_matrix(17,1).*exp(-efficiency_matrix(17,2).*x.^efficiency_matrix(17,3)) .* (1-exp(-efficiency_matrix(17,4).*x.^efficiency_matrix(17,5)))  );
        end     
    end
    
    % Turn me off!!!  Just for testing Hannah output
    %effcal=1;
    
    % Number of columns - testing shelf regex
%     num_cols = 10;
    num_cols = 9;
    
    
    extract_energies = C3{1,1};
    dim = length(C3{1,1}); 
    
    % Make matrix to hold data for this spectrum
    sto_mat = zeros(dim,num_cols);
    
    
    % Platform: 1 Type:
    % pattern = 'Platform:\s+(\d+)';
    % data = regexp(indata, pattern, 'tokens');
    % data = cell2mat(cellfun(@(x) str2double(x{:}), data, 'UniformOutput', false));
    
    % Extract header
    mass = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens');
    mass = cell2mat(cellfun(@(x) str2double(x{:}), mass, 'UniformOutput', false));
    mass_str = num2str(mass);
    foil_id = str2num(mass_str(1:end-2));
    
    ct = regexp(raw_str, 'Datetime:\s+([-\w: ]+)', 'tokens');
%     ct{1,1}{1,1}
%     ct = cell2mat(cellfun(@(x) str2double(x{:}), ct, 'UniformOutput', false))
    ct = juliandate(datetime(ct{1,1}{1,1}));
    
    cl = regexp(raw_str, 'Live:\s+(\d+)', 'tokens');
    cl = cell2mat(cellfun(@(x) str2double(x{:}), cl, 'UniformOutput', false));
    
    % Append data from header to columns
    sto_mat(:,4) = mass;
    sto_mat(:,5) = ct;
    sto_mat(:,6) = cl;
    
    % t12 = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens');
    % t12 = cell2mat(cellfun(@(x) str2double(x{:}), t12, 'UniformOutput', false))
    %
    % Igamma = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens');
    % Igamma = cell2mat(cellfun(@(x) str2double(x{:}), Igamma, 'UniformOutput', false))

    
    % Pull out gamma-ray energies and net counts
%                     energy        net counts (efficiency and attenuation corrected)                                      %error
    sto_mat(:,1:3) = [C3{1,1}      double(C3{1,6})./(effcal(C3{1,1}).*exp(-ppval(pp2,C3{1,1}).* 0.5e-3 .* rhodrs(foil_id)) )       C3{1,7}];
%     sto_mat(:,1:3) = [C3{1,1}      double(C3{1,6})./(effcal.*         exp(-ppval(pp2,C3{1,1}).* 0.5e-3 .* rhodrs(foil_id)) )       C3{1,7}];
%     sto_mat(:,10) = shelf;
    
    
    % Replace energy values with closest from key energies, otherwise delete?
    % Indices of rows to keep
    ind = zeros(2,dim);
    
    for j=1:dim
        [Ydiff, idx] = min(abs(key_energies - sto_mat(j,1)));
        if Ydiff <= 1.0
            ind(:,j) = [j;idx];
        end
    end
    
    % Remove null values
    % ind = ind(ind ~= 0);
    ind = ind(:,ind(1,:) ~= 0);
    
    % Extract rows within 1 keV or key energies
    sto_mat = sto_mat(ind(1,:), :);
    % Replace with key energies
    sto_mat(:,1) = key_energies(ind(2,:));
    % append the t12 and Igamma data from a key
    sto_mat(:,7:9) = glines(ind(2,:),:);
    
    
    
    % Merge matrix into data_mat
    
    % [~,ii] = ismember(A(:,2),B(:,1));
    % [~,ii] = ismember(sto_mat(:,1),key_energies(:,1));
    % out = [key_energies, sto_mat(ii,2)];
    
    % data_mat(ind(2,:),(2:7)+7.*(i-1)+1) = sto_mat;
    
    % Need to append the t12 and Igamma data from a key ?
    % i=1;
%     data_mat(ind(2,:),(1:8)+8.*(i-1)+1) = sto_mat(:,2:length(sto_mat));
    data_mat(ind(2,:),(1:8)+8.*(i-1)+1) = sto_mat(:,2:9);

    % data_mat(ind(2,:),((1:5)+5.*(i-1)+2):((1:5)+5.*(i-1)+4)) = sto_mat(:,2:length(sto_mat))
    
end

unique(shelves)



fclose('all');

