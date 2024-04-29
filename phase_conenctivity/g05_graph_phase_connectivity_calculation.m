function g05_graph_phase_connectivity_calculation(paths,graph_ps)

input_root_path = paths.tf_decomposed_data_4phase_connectivity;
output_root_path = paths.phase_connectivity_results_raw;
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
    all2all_phrase = connectivity_calculation(tf_phrase_data);
    
    %% sentence data ...
    tmp_sentence_data_name = [tmp_subj,'_tf_graph_sentence.mat'];
    load(fullfile(tmp_data_path,tmp_sentence_data_name)); % tf_sentence_data ...
    fprintf('%s(sentence) connectivity calculation ... \n\n', tmp_subj);
    all2all_sentence = connectivity_calculation(tf_sentence_data);
    
    %% concatenate all2all connectivity data ...
    fprintf('%s binarizing connectivity results ... \n\n ', tmp_subj);
    [conn_degree_phrase_raw,conn_degree_sentence_raw,...
        connmat_phrase_raw,connmat_sentence_raw,...
        connmat_phrase_binarized,connmat_sentence_binarized]...
        = graph_thresholding(all2all_phrase,all2all_sentence,n_std_above);
    
    %% saving results ...
    tmp_conn_results_saving_path = fullfile(output_root_path,tmp_subj);
    if ~exist(tmp_conn_results_saving_path,'dir')
        mkdir(tmp_conn_results_saving_path);
    end
    
    save(fullfile(tmp_conn_results_saving_path,[tmp_subj,'_phrase_connectivity_results.mat']),...
        'conn_degree_phrase_raw',...
        'connmat_phrase_raw',...
        'connmat_phrase_binarized',...
        'freq_range','times');
    save(fullfile(tmp_conn_results_saving_path,[tmp_subj,'_sentence_connectivity_results.mat']),...
        'conn_degree_sentence_raw',...
        'connmat_sentence_raw',...
        'connmat_sentence_binarized',...
        'freq_range','times');
    
end

function tf_all2all = connectivity_calculation(tf_data)
%% nchan*nchan*frex*times ...
tf_all2all = zeros(size(tf_data,1),size(tf_data,1),size(tf_data,2),size(tf_data,3));

%% now that we have all the data, compute all-to-all connectivity
for chani=1:size(tf_data,1)
    for chanj=chani+1:size(tf_data,1)
        
        % compute connectivity (cross spectral density) ...
        xsd = squeeze(tf_data(chani,:,:,:).*conj(tf_data(chanj,:,:,:)));
        % connectivity matrix  ...
        tf_all2all(chani,chanj,:,:) = abs(mean(exp(1i*angle(xsd)),3)); % ispc
        %         tf_all2all(chani,chanj,:,:) = abs(mean(sign(imag(xsd)),3)); % pli
        
    end
end

function [conn_degree_phrase_raw,conn_degree_sentence_raw,...
    connmat_phrase_raw,connmat_sentence_raw,...
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
% % phrase_base_avg = repmat(mean(mean(conn_degree_phrase_raw(:,:,baseidx(1):baseidx(2)),3),2),[1 size(conn_degree_phrase_raw,2) size(conn_degree_phrase_raw,3)]);
% phrase_base_avg(phrase_base_avg == 0) = 1e-2;
% % conn_degree_phrase_base = conn_degree_phrase_raw - phrase_base_avg;
% conn_degree_phrase_base_percent = 100*(conn_degree_phrase_raw - phrase_base_avg)./phrase_base_avg;
% 
% %% subtract baseline for sentence ...
% sentence_base_avg = repmat(mean(conn_degree_sentence_raw(:,:,baseidx(1):baseidx(2)),3),[1 1 size(conn_degree_sentence_raw,3)]);
% % sentence_base_avg = repmat(mean(mean(conn_degree_sentence_raw(:,:,baseidx(1):baseidx(2)),3),2),[1 size(conn_degree_sentence_raw,2) size(conn_degree_sentence_raw,3)]);
% sentence_base_avg(sentence_base_avg == 0) = 1e-2;
% % conn_degree_sentence_base = conn_degree_sentence_raw - sentence_base_avg;
% conn_degree_sentence_base_percent = 100*(conn_degree_sentence_raw - sentence_base_avg)./sentence_base_avg;
