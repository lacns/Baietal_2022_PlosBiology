function p05_pac_z_calculation_and_z_transformation(paths,pac_ps)
input_root_path = paths.power_phase_time_series; 
output_root_path = paths.pacz;

%% listing all subjects ...
all_subjects = dir(input_root_path);
all_subjects(1:2)=[];
n_iter = pac_ps.n_iter;

%% calculation ...
for subj_i = 6:length(all_subjects) %##############
    tmp_subj = all_subjects(subj_i).name;
    tmp_pac_path = fullfile(input_root_path,tmp_subj);
    %% phrase condition  ...
    tmp_pac_phrase_name = [tmp_subj,'_pac_ts_phrase.mat'];
    load(fullfile(tmp_pac_path,tmp_pac_phrase_name));
    tmp_power_data_phrase = phrase_ts.power;
    tmp_phase_data_phrase = phrase_ts.phase;
    clearvars phrase_ts;
    pacz_phrase = [];
    for chan_i = 1:length(tmp_power_data_phrase)
        tmp_chan_power_ts_phrase = tmp_power_data_phrase(chan_i,:);
        tmp_chan_phase_ts_phrase = tmp_phase_data_phrase(chan_i,:);
        tmp_pacz_phrase = ...
            pacZ_calculation(tmp_chan_power_ts_phrase,tmp_chan_phase_ts_phrase,n_iter,...
            tmp_subj,'phrase',chan_i);
        pacz_phrase(:,:,chan_i) = tmp_pacz_phrase;
    end
    
    %% sentence condition ...
    tmp_pac_sentence_name = [tmp_subj,'_pac_ts_sentence.mat'];
    load(fullfile(tmp_pac_path,tmp_pac_sentence_name));
    tmp_power_data_sentence = sentence_ts.power;
    tmp_phase_data_sentence = sentence_ts.phase;
    clearvars sentence_ts;
    pacz_sentence = [];
    for chan_i = 1:length(tmp_power_data_sentence)
        tmp_chan_power_ts_sentence = tmp_power_data_sentence(chan_i,:);
        tmp_chan_phase_ts_sentence = tmp_phase_data_sentence(chan_i,:);
        tmp_pacz_sentence = ...
            pacZ_calculation(tmp_chan_power_ts_sentence,tmp_chan_phase_ts_sentence,n_iter,...
            tmp_subj,'sentence',chan_i);
        pacz_sentence(:,:,chan_i) = tmp_pacz_sentence;
    end
    
    %% saving results ...
    tmp_pacz_saving_path = fullfile(output_root_path,tmp_subj);
    if ~exist(tmp_pacz_saving_path,'dir')
        mkdir(tmp_pacz_saving_path);
    end
    save(fullfile(tmp_pacz_saving_path,[tmp_subj,'_pacZ_phrase.mat']),'pacz_phrase');
    save(fullfile(tmp_pacz_saving_path,[tmp_subj,'_pacZ_sentence.mat']),'pacz_sentence');
    clc;
end

function pacz = pacZ_calculation(power_ts,phase_ts,n_iter,tmp_subj,cond_label,chan_i)
rng('shuffle');
n_bins = length(power_ts)*length(phase_ts);
idx=1;
for power_i = 1:length(power_ts)
    tmp_freq_power_ts = power_ts{power_i};
    
    for phase_i = 1:length(phase_ts)
        fprintf(['%s, condition %s, channel %d \n', 'pac z-transformation (power band:%d , phase band: %d), %d of %d ...\n\n'],...
            tmp_subj,cond_label,chan_i,power_i,phase_i,idx,n_bins);
        tmp_freq_phase_ts = phase_ts{phase_i};
        %% compute observed PAC
        obsPAC = abs(mean( reshape(tmp_freq_power_ts,1,[]).*exp(1i*reshape(tmp_freq_phase_ts,1,[])) ));
        
        %% PAC z-transfprmation by permutaions ...
        n_trials = size(tmp_freq_phase_ts,2);
        data_length = size(tmp_freq_phase_ts,1);
        %         n_iter = 500; % n iterations ...
        for i_iter=1:n_iter
            for trial_i=1:n_trials
                cut_point = randperm(data_length,1);
                tmp_freq_power_ts_null(:,trial_i) = tmp_freq_power_ts([cut_point:end, 1:cut_point-1],trial_i);
            end
            permutedPAC(i_iter) = abs(mean( reshape(tmp_freq_power_ts_null,1,[]).*exp(1i*reshape(tmp_freq_phase_ts,1,[])) ));
        end
        pacz(power_i,phase_i) = (obsPAC-mean(permutedPAC))/std(permutedPAC);
        idx = idx+1;
    end
end