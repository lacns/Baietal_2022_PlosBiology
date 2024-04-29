function a03_auditory_stage_tf_decomposition(paths,tf_ps)
input_root_path = paths.condition_seperated_data;
output_root_path = paths.tf_results_raw;

all_subjects = dir(paths.condition_seperated_data);
all_subjects(1:2) = [];

%%
for subj_i=1:length(all_subjects)
    tmp_subj = all_subjects(subj_i).name;
    tmp_subj_data_path = fullfile(input_root_path,all_subjects(subj_i).name);
    
    %% phrase listening (auditory stage)...
    tmp_subj_phrase_data_name = ['c00_',tmp_subj,'_clean_csd_auditory_phrases.set'];
    tmp_eeg_phrase = pop_loadset('filename',tmp_subj_phrase_data_name,'filepath',tmp_subj_data_path);
    tmp_subj_phrase_tf_results = wavelet_convolution(tmp_eeg_phrase,tf_ps,tmp_subj,'phrase');
    
    %% sentence listening (auditory stage)...
    tmp_subj_sentence_data_name = ['c00_',tmp_subj,'_clean_csd_auditory_sentences.set'];
    tmp_eeg_sentence = pop_loadset('filename',tmp_subj_sentence_data_name,'filepath',tmp_subj_data_path);
    tmp_subj_sentence_tf_results = wavelet_convolution(tmp_eeg_sentence,tf_ps,tmp_subj,'sentence');
    
    %% saving data ...
    tmp_subj_tf_saving_path = fullfile(output_root_path,tmp_subj);
    if ~exist(tmp_subj_tf_saving_path,'dir')
        mkdir(tmp_subj_tf_saving_path);
    end
    fprintf('saving %s tf decomposition results ... \n\n', tmp_subj);
    save(fullfile(tmp_subj_tf_saving_path,[tmp_subj,'_phrase_tf_results_auditory_stage.mat']),'tmp_subj_phrase_tf_results');
    save(fullfile(tmp_subj_tf_saving_path,[tmp_subj,'_sentence_tf_results_auditory_stage.mat']),'tmp_subj_sentence_tf_results');
    close all; clc;
end

function tf_results = wavelet_convolution(eeg,tf_ps,tmp_subj,audio_type)
%% wavelet parameters
n_frex = tf_ps.n_frex;
min_freq =  tf_ps.freq_min;
max_freq = tf_ps.freq_max;
range_cycles = tf_ps.cycle_range;
% baseline = tf_ps.baseline;


frex = linspace(min_freq,max_freq,n_frex);
times = -1:1/eeg.srate:1;
half_wave = (length(times)-1)/2;


%% FFT parameters
nKern = length(times);
nData = eeg.pnts*eeg.trials;
nConv = nKern+nData-1;

%% convert baseline time into indices
% [~,baseidx(1)] = min(abs(eeg.times-baseline(1)));
% [~,baseidx(2)] = min(abs(eeg.times-baseline(2)));

%% set a few different wavelet widths (number of wavelet cycles)
% range_cycles = [3 30]; % [2.5  50];

%% other wavelet parameters
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),n_frex);
% nCycles = linspace(range_cycles(1),range_cycles(2),num_frex);

for chan_i = 1:eeg.nbchan
    
    %% FFT of data (doesn't change on frequency iteration)
    channel2use = eeg.chanlocs(chan_i).labels;
    dataX = fft(reshape(eeg.data(strcmpi(channel2use,{eeg.chanlocs.labels}),:,:),1,eeg.pnts*eeg.trials),nConv);
    
    %% initialize output time-frequency data
    tf_power = zeros(length(frex),eeg.pnts);
    tf_itpc = zeros(n_frex,eeg.pnts);
    
    %% frequency loop ...
    for fi=1:length(frex)
        
        %% create wavelet and get its FFT
        s = nCycles(fi)/(2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*times) .* exp(-times.^2./(2*s^2));
        cmwX = fft(cmw,nConv);
        
        %% run convolution
        as = ifft(cmwX.*dataX,nConv);
        as = as(half_wave+1:end-half_wave);
        as = reshape(as,eeg.pnts,eeg.trials);
        
        %% put power and itpc data into big matrix
        tf_power(fi,:) = mean(abs(as).^2,2);
        tf_itpc(fi,:) = abs(mean(exp(1i*angle(as)),2));
        
    end
    fprintf('%s TF decomposition (auditory: %s): %s (number: %d)... \n\n',tmp_subj,audio_type,channel2use,chan_i);
    
    %% db conversion
%     tf_power = 10*log10( bsxfun(@rdivide, tf_power, mean(tf_power(:,baseidx(1):baseidx(2)),2)) );
    %     tf_itpc = 10*log10( bsxfun(@rdivide, tf_itpc, mean(tf_itpc(:,baseidx(1):baseidx(2)),2)) );
    tf_results.tf_power(chan_i,:,:) = tf_power;
    tf_results.tf_itpc(chan_i,:,:) = tf_itpc;
    
end
tf_results.times = eeg.times;
tf_results.frex = frex;
tf_results.cycle = nCycles;
tf_results.chan_locs = eeg.chanlocs;
