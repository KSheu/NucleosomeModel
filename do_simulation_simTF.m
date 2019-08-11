
% TNF SIMULATION
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

global START_TIME END_TIME 


START_TIME =0;
END_TIME = 480;
% % END_TIME = 900;
% Transcription factor vector


% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM',  'nfkb_curves_LPS100ng','nfkb_curves_pic50ug'};
% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM'};
% names = {'nfkb_curves_TNF10ng','nfkb_curves_TNF_ikbamm'};
% names = {'nfkb_oscillatory','nfkb_nonoscillatory', 'nfkb_oscillatory_hiamp', 'nfkb_nonoscillatory_hiamp'};
names = {'nfkb_oscillatory','nfkb_nonoscillatory'};
% names = {'nfkb_oscillatory_2xtotalactivity','nfkb_persistent_2xtotalactivity'}; %use END_TIME=900 if this TF sim
for j = 1:length(names)
data_name = char(names(j));
% data_name = 'nfkb_curves_TNF';
data = load(strcat('F://enhancer_dynamics/nfkb_trajectories/simTFs/',data_name,'.mat'));

data = cell2struct(struct2cell(data), {'nfkb_curves'});
% data = (data.nfkb_curves)*30;
data = (data.nfkb_curves)*1;



% time = START_TIME:END_TIME;
% tf = data_smooth(5,1:96); %cut to 8hrs
% plot(sim.time*5, interp1( 1:length(tf), tf, START_TIME:END_TIME,'nearest'));
% xlabel('time (min)')
% ylabel('[TF]')	
% ylim([0 inf]);
% f = fit( transpose(time(:,1:96)*5), transpose(tf), 'linearinterp');
% plot( f, transpose(time(:,1:96)*5), transpose(tf) )

% Starting Conditions
initvalues = zeros(15,1);
initvalues(1,1) = 1;    %E0
initvalues(2,1) = 0;    %E1

initvalues(3,1) = 0;   %E2
initvalues(4,1) = 0;   %E3

initvalues(5,1) = 0;   %E4
initvalues(6,1) = 0;   %E5
initvalues(7,1) = 0;   %E6
initvalues(8,1) = 0;   %E7
initvalues(9,1) = 0;   %E8
initvalues(10,1) = 0;   %E9
initvalues(11,1) = 0;   %E10
initvalues(12,1) = 0;   %E11
initvalues(13,1) = 0;   %E12
initvalues(14,1) = 0;   %E13
initvalues(15,1) = 0;   %E14



tf = transpose(data); %cut to 8hrs
time = linspace(START_TIME, END_TIME, length(tf));
[tsim1, results1] = ode15s(@(t,y) chromatinOde(t, y, time,{}, tf),[START_TIME END_TIME],initvalues);

%   [t_sim, x_sim] = ode15s(@(t,y)'chromatinOde', v.SIM_TIME,starting_vals,ode_opt,v);

output = transpose(interp1(tsim1,results1,START_TIME:END_TIME, 'linear'));

subplot(1,length(names),j);
plot(output(15,:));
ylim([0 1]);
title(char(names(j)));
colorbar
end



% plot(output(15,:));

