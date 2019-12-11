clear all; clc;

% Load parsed matrices from ekeberg_spec_anal.m
load ni_mat.mat;
load fe_mat.mat;
load cu_mat.mat;
% load ir_mat.mat;


mat = ni_mat;
bad_files =  [],

for i=128:100:1028
     [~,col] = find((mat == i)  )
     bad_files = [bad_files ; col]
end
for i=129:100:1029
     [~,col] = find((mat == i)  )
     bad_files = [bad_files ; col]
end
% for i=126:100:1026
%      [~,col] = find((mat == i)  )
%      bad_files = [bad_files ; col]
% end
for i=177:100:1077
     [~,col] = find((mat == i)  )
     bad_files = [bad_files ; col]
end

bad_files = unique(bad_files);

bad_cts = mat(:,bad_files+1);
bad_masses = mat(:,bad_files);


% EoB Time
delta_ts = (bad_cts - juliandate(datetime('26-Feb-2019 00:32:00')) ) .* 24;
