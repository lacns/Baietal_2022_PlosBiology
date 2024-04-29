function p05_time_sliding_pacz_coupling_calculation(paths,pac_ps)

%% listing all subjects ...
all_subjects = dir(paths.eeg_after_ht_path);
all_subjects(1:2) = [];

for subj_i = 1:length(all_subjects)
    tmp_subj = all_subjects(subj_i).name;
    tmp_data_path = fullfile(paths.eeg_after_ht_path,tmp_subj);
    
    %% phrase data ...
    tmp_phrase_data_name = [tmp_subj,'_ts_phrase.mat'];
    load(fullfile(tmp_data_path,tmp_phrase_data_name));
    eeg_phrase = pacz_calculation(eeg_phrase,pac_ps,tmp_subj,'phrase');
    
    %% sentence data ...
    tmp_sentence_data_name = [tmp_subj,'_ts_sentence.mat'];
    load(fullfile(tmp_data_path,tmp_sentence_data_name));
    eeg_sentence = pacz_calculation(eeg_sentence,pac_ps,tmp_subj,'sentence');
    
    %% saving data ...
    pacz_results_saving_path = fullfile(paths.pacz_results_path,tmp_subj);
    if ~exist(pacz_results_saving_path,'dir')
        mkdir(pacz_results_saving_path);
    end
    
    fprintf('saving %s data ...\n\n',tmp_subj);
    save(fullfile(pacz_results_saving_path,[tmp_subj,'_pacz_results_phrase.mat']),'eeg_phrase');
    save(fullfile(pacz_results_saving_path,[tmp_subj,'_pacz_results_sentence.mat']),'eeg_sentence');
    
end

function eeg = pacz_calculation(eeg,pac_ps,tmp_subj,data_type)
win_centers = dsearchn(eeg.times',pac_ps.ts_center');
win_length= (dsearchn(eeg.times',pac_ps.ts_length) - dsearchn(eeg.times',0))-1; % win length - center ...
win_length_left = ceil(win_length/2);
win_length_right = win_length - win_length_left;
n_iter = pac_ps.n_iter;
for time_i = 1:length(win_centers)
    tmp_center = win_centers(time_i);
    tmp_win = [tmp_center - win_length_left,tmp_center+win_length_right];
%     tmp_win = [tmp_center, tmp_center+win_length];
    fprintf('%s pacz calculation (%s) at time %d ms, (time idx from: %d of %d) ...\n\n',...
        tmp_subj,data_type,round(pac_ps.ts_center(time_i)),time_i,length(win_centers));
    tmp_phase_ts = eeg.phase_ts(:,tmp_win(1):tmp_win(2),:);
    tmp_power_ts = eeg.power_ts(:,tmp_win(1):tmp_win(2),:);
    tmp_pacz = calculation_and_permutation(tmp_phase_ts,tmp_power_ts,n_iter);
    eeg.pacz(:,time_i) = tmp_pacz;
    eeg.pacz_time(1,time_i) = round(pac_ps.ts_center(time_i));
end

function pacz = calculation_and_permutation(phase_ts,power_ts,n_iter)
rng('shuffle');
for chan_i = 1:size(phase_ts,1)
    tmp_chan_phase_ts = squeeze(phase_ts(chan_i,:,:));
    tmp_chan_power_ts = squeeze(power_ts(chan_i,:,:));
    obsPAC = abs(mean((reshape(tmp_chan_power_ts,1,[]).*exp(1i*reshape(tmp_chan_phase_ts,1,[])))));
    
    %% permutations ...
    n_trials = size(tmp_chan_power_ts,2); % for each trial shift the same amount ...
    data_length = size(tmp_chan_power_ts,1);
    for i_iter=1:n_iter
        for trial_i=1:n_trials
            cut_point = randperm(data_length,1);
            % only shift power time series ...
            shifted_tmp_power_ts(:,trial_i) = tmp_chan_power_ts([cut_point:end, 1:cut_point-1],trial_i);
        end
        permutedPAC(i_iter) = abs(mean( reshape(shifted_tmp_power_ts,1,[]).*exp(1i*reshape(tmp_chan_phase_ts,1,[])) ));
    end
    pacz(chan_i,1) = (obsPAC-mean(permutedPAC))/std(permutedPAC);
end


