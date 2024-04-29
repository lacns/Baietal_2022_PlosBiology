function ls08_strf_model_fitting(paths,strf_ps,lambda_info)
%% adding path ...
addpath(paths.mtrf_toolbox);
input_root_path = paths.strf_sets_seperated_data;
output_root_path = paths.strf_model_results;

%% training ...
all_subjects = dir(fullfile(input_root_path,'*strf_cv_dataset*.mat'));
for subj_i=1:length(all_subjects)
    
    %%
    load(fullfile(input_root_path,all_subjects(subj_i).name)); % strf_cv ...
    
    %% training data for tmp subject ...
    tmp_training_data = [strf_cv.training_data_phrase,strf_cv.training_data_sentence];
    
    %% best lambda for tmp subject  ...
    tmp_subj = all_subjects(subj_i).name(1:8);
    tmp_subj_idx = strcmpi({lambda_info.subj_id},tmp_subj);
    tmp_lambda_info = lambda_info(tmp_subj_idx);
    
    %% bands bins ...
    unique_bands = unique({tmp_lambda_info.band_name});
    
    %% model fitting ...
    for trial_i = 1:length(tmp_training_data)
        tmp_stim = tmp_training_data(trial_i).spectrogram; % tmp_stim ...
        tmp_audio_type = tmp_training_data(trial_i).audio_type;
        
        for band_i = 1:length(unique_bands)
            tmp_band_type = unique_bands{band_i};
            if isempty(tmp_band_type)
                tmp_resp_name = 'eeg_data';
                tmp_model_name = 'model_one_band';
            else
                tmp_resp_name = ['eeg_data_',tmp_band_type]; % eeg data ...
                tmp_model_name = ['model_',unique_bands{band_i}];
            end
            
            eval(['tmp_resp = tmp_training_data(',num2str(trial_i),').',tmp_resp_name,';']); % tmp_resp ...
            tmp_lambda_idx = strcmpi({tmp_lambda_info.strf_type},tmp_audio_type) & ...
                strcmpi({tmp_lambda_info.band_name},tmp_band_type);
            tmp_lambda = tmp_lambda_info(tmp_lambda_idx).best_lambda;
            
            % direction = 1 ... 
            tmp_model = mTRFtrain(tmp_stim,tmp_resp,...
                strf_ps.fs,1,strf_ps.training_tmin,strf_ps.training_tmax,tmp_lambda,'method','ridge',...
                'zeropad',0);
            eval(['tmp_training_data(',num2str(trial_i),').',tmp_model_name,' = tmp_model;']);
            fprintf('%s (%s, %s) forward model fitting trial %d of %d... \n\n',...
                tmp_subj,tmp_audio_type,tmp_model_name,trial_i,length(tmp_training_data));
        end
    end
    
    %% saving fitted model for each subject ...
    strf_model_results.training_data = tmp_training_data;
    strf_model_results.testing_data = [strf_cv.testing_data_phrase,strf_cv.testing_data_sentence];
    save(fullfile(output_root_path,[tmp_subj,'_strf_fitting_results.mat']),'strf_model_results','chan_locs');
end