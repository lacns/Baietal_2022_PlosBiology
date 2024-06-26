function [paths,dora_ps,idx] = d01_path_and_parameters_setting(root_path)
%% path setting ... 
method_specific_path = '\step_03_eeg_data_analysis\e06_eeg_results';
paths.result_path = fullfile(root_path,method_specific_path,'r08_dora_simulation');
making_folder(paths.result_path);

%% task index ... 
idx.define_units = 0;%# ...
idx.generatng_stimuli = 0;%# ... 
idx.simulation = 0;%# ... 
idx.freq_decomposition = 0;%# ...
idx.statistics = 0;%# ...
idx.plot = 1;%# ...

%% parameters setting ... 
dora_ps.n_trials = 100; % 100 trials ... 
dora_ps.syllable_win = 250; % 250 ms for each syllable ...
dora_ps.n_syllable_per_sec = 4; % 4 sylables per second ...
dora_ps.trial_length = 12; % whole trial length 12-second ... 
dora_ps.fs = 1000; % in ms ...
dora_ps.syllable_range = [140, 230]; % syllable length for each syllable in ms >>> corresponds to 4 to 7 Hz ... 
dora_ps.buffer = 3; % in point >>> syllable started after 3 ms and ended 3 ms before the end ... 
dora_ps.roi_units = 'po_unit';
dora_ps.n_smooth = 35; % signal construction ... 
dora_ps.n_binding = [30,250,350];
dora_ps.freq_limit = 13.5; % 6.5 Hz below ...
dora_ps.n_itpc = 100; % times for calculating itpc ... 
dora_ps.n_per_calculation = 50;
dora_ps.coupling4phrase = {'f','c'};
dora_ps.coupling4sentence = {'g','c'};
dora_ps.n_c = 100;
dora_ps.n_c_per_calculation = 30;
dora_ps.order_filtering = [50,35,20];


%% creating folders ... 
function making_folder(folder_name)
if ~exist(folder_name,'dir')
mkdir(folder_name);
end

