function g05_graph_power_connectivity_calculation(paths,graph_ps)

input_root_path = paths.tf_decomposed_data_4power_connectivity;
output_root_path = paths.power_connectivity_results_raw;
all_subjects = dir(input_root_path);
all_subjects(1:2) = [];

freq_range = graph_ps.frex;
times = graph_ps.times2save;
n_std_above = graph_ps.n_std_above;
% baseidx = graph_ps.baseidx;
%%
for subj_i = 1:length(all_subjects)
    tmp_subj = all_subjects(subj_i).name;
    tmp_data_path = fullfile(input_root_path,tmp_subj);
    
    %% phrase data ...
    tmp_phrase_data_name = [tmp_subj,'_tf_graph_phrase.mat'];
    load(fullfile(tmp_data_path,tmp_phrase_data_name)); % tf_phrase_data ...
    fprintf('%s(phrase) connectivity calculation ... \n\n', tmp_subj);
    all2all_phrase = connectivity_calculation(tf_phrase_data,graph_ps);
    %     all2all_phrase = atanh(all2all_phrase); % Fisher-Z transformation ...
    
    %% sentence data ...
    tmp_sentence_data_name = [tmp_subj,'_tf_graph_sentence.mat'];
    load(fullfile(tmp_data_path,tmp_sentence_data_name)); % tf_sentence_data ...
    fprintf('%s(sentence) connectivity calculation ... \n\n', tmp_subj);
    all2all_sentence = connectivity_calculation(tf_sentence_data,graph_ps);
    %     all2all_sentence = atanh(all2all_sentence); % Fisher-Z transformation ...
    
    %% concatenate all2all connectivity data ...
    fprintf('%s binarizing connectivity results ... \n\n ', tmp_subj);
    [conn_degree_phrase_raw,conn_degree_sentence_raw,...
        connmat_phrase_binarized,connmat_sentence_binarized]...
        = graph_thresholding(all2all_phrase,all2all_sentence,n_std_above);
    
    %% saving results ...
    tmp_conn_results_saving_path = fullfile(output_root_path,tmp_subj);
    if ~exist(tmp_conn_results_saving_path,'dir')
        mkdir(tmp_conn_results_saving_path);
    end
    
    save(fullfile(tmp_conn_results_saving_path,[tmp_subj,'_phrase_power_connectivity_results.mat']),...
        'conn_degree_phrase_raw',...
        'connmat_phrase_binarized',...
        'freq_range','times');
    save(fullfile(tmp_conn_results_saving_path,[tmp_subj,'_sentence_power_connectivity_results.mat']),...
        'conn_degree_sentence_raw',...
        'connmat_sentence_binarized',...
        'freq_range','times');
    
end

function power_corr_all2all = connectivity_calculation(tf_data,graph_ps)
%% nchan*nchan*frex*times ...
power_corr_all2all = zeros(size(tf_data,1),size(tf_data,1),size(tf_data,2),size(tf_data,3));
tf_power = abs(tf_data).^2;
baseline_idx = graph_ps.baseidx;
% avg_baseline_over_trials = squeeze(mean(tf_power,4));
% avg_baseline = squeeze(mean(avg_baseline_over_trials(:,:,graph_ps.baseidx(1):graph_ps.baseidx(2)),3));
% tf_power_base_removed = 10*log10(tf_power./repmat(avg_baseline,1,1,size(tf_power,3),size(tf_power,4)));

%% now that we have all the data, compute all-to-all connectivity
tic;
% for chan_i=1:size(tf_power,1)
%     for chan_j=chan_i+1:size(tf_power,1)
%         for freq_i = 1:size(tf_power,2)
%             for time_i  = 1:size(tf_power,3)
%                 tmp_signal1 = squeeze(tf_power(chan_i,freq_i,time_i,:));
%                 tmp_signal2 = squeeze(tf_power(chan_j,freq_i,time_i,:));
%                 power_corr_all2all(chan_i,chan_j,freq_i,time_i) = corr(tmp_signal1,tmp_signal2,'type','s');
%             end
%         end
%     end
% end
% n_chan = size(tf_power,1);
n_frex = size(tf_power,2);
n_times = size(tf_power,3);
n_trials = size(tf_power,4);

for chan_i=1:size(tf_power,1)
    tmp_chan_i_data = reshape(permute(squeeze(tf_power(chan_i,:,:,:)),[3,1,2]),n_trials,[]);    
    for chan_j=chan_i+1:size(tf_power,1)
        tmp_chan_j_data = reshape(permute(squeeze(tf_power(chan_j,:,:,:)),[3,1,2]),n_trials,[]);
        tmp_corr = arrayfun(@(idx) corr(tmp_chan_i_data(:,idx), tmp_chan_j_data(:,idx),'type','Spearman'), 1:n_frex*n_times, 'Uni', 1);
        tmp_corr = reshape(tmp_corr,n_frex,n_times);
        power_corr_all2all(chan_i,chan_j,:,:) = tmp_corr;
    end
end
toc;
%%
function [conn_degree_phrase_raw,conn_degree_sentence_raw,...
    connmat_phrase_binarized,connmat_sentence_binarized]= ...
    graph_thresholding(tf_all2all_phrase,tf_all2all_sentence,n_std_above)
n_frex = size(tf_all2all_phrase,3);
n_times = size(tf_all2all_phrase,4);
conn_all2all = cat(5,tf_all2all_phrase,tf_all2all_sentence);

%% threshold the connectivity matrix (separate threshold for each frequency)
for fi=1:n_frex
    
    %% define threshold ...
    tempsynch = nonzeros(conn_all2all(:,:,fi,:,:)); % threshold based on all time bins ...
    thresh(fi) = median(tempsynch) + n_std_above*std(tempsynch);
    
    % isolate, threshold, binarize
    for ti=1:n_times
        %% for phrase ...
        temp_phrase = squeeze(tf_all2all_phrase(:,:,fi,ti));
        temp_phrase = temp_phrase + triu(temp_phrase)';      % make symmetric matrix
        connmat_phrase_raw(:,:,fi,ti) = temp_phrase;
        
        temp_phrase = temp_phrase > thresh(fi);         % threshold and binarize
        connmat_phrase_binarized(:,:,fi,ti) = temp_phrase;
        
        % number of connections in each time&frequency point ...
        conn_degree_phrase_raw(:,fi,ti) = sum(temp_phrase); % compute degree (sum of suprathreshold connections)
        
        %% for sentence ...
        temp_sentence = squeeze(tf_all2all_sentence(:,:,fi,ti));
        temp_sentence = temp_sentence + triu(temp_sentence)';      % make symmetric matrix
        connmat_sentence_raw(:,:,fi,ti) = temp_sentence;
        
        temp_sentence = temp_sentence > thresh(fi);         % threshold and binarize
        connmat_sentence_binarized(:,:,fi,ti) = temp_sentence;
        
        % number of connections in each time&frequency point ...
        conn_degree_sentence_raw(:,fi,ti) = sum(temp_sentence);
        
    end
end

% %%  subtract baseline for phrase ...
% phrase_base_avg = repmat(mean(conn_degree_phrase_raw(:,:,baseidx(1):baseidx(2)),3),[1 1 size(conn_degree_phrase_raw,3)]);
% phrase_base_avg(phrase_base_avg == 0) = 1;
% conn_degree_phrase_base_percent = 100*(conn_degree_phrase_raw - phrase_base_avg)./phrase_base_avg;
% 
% %% subtract baseline for sentence ...
% sentence_base_avg = repmat(mean(conn_degree_sentence_raw(:,:,baseidx(1):baseidx(2)),3),[1 1 size(conn_degree_sentence_raw,3)]);
% sentence_base_avg(sentence_base_avg == 0) = 1;
% conn_degree_sentence_base_percent = 100*(conn_degree_sentence_raw - sentence_base_avg)./sentence_base_avg;
% 
