names = {'TNF10ng_754','LPS100ng_756', 'P3CSK4100ng_547','CpG100nM_610', 'PIC50ug_566',};
% names = {'WT_10ngTNF', 'ikbamut_10ngTNF'};
for j = 1:length(names)
    
data_name = char(names(j));
data = load(strcat('F://enhancer_dynamics/nfkb_trajectories_08142019/nfkb_dynamics_',data_name,'.mat'));

data = cell2struct(struct2cell(data), {'trajectories'});
data = table2array(data.trajectories);
data_smooth = smoothrows(data);
data_smooth(any(isnan(data_smooth), 2), :) = []; %remove NaN rows

datazero = vec2mat(data_smooth(:, 1), 1);
subtract = repmat(datazero, 1, size(data,2));
data_smooth = data_smooth - subtract; %subtract the first column to normalize
data_smooth(data_smooth<0) = 0; %take neg. values to be 0

save(strcat('F://enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_',data_name,'_smoothed.mat'), 'data');

end
