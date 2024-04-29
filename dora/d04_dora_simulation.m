function d04_dora_simulation(paths, dora_ps)
input_path = paths.result_path;
output_path = input_path;

%% loading units' definition ...
load(fullfile(input_path,'units_definition'));

%% loading simulated trials ...
fprintf('loading simulated input ... \n');
load(fullfile(input_path,'simulated_stimuli'));
n_str = fprintf('phrase simulations (%3d of %3d) ... ',0,0);

%% 1st >>> phrase simulation ...
for i=1:dora_ps.n_trials
    %% syllable input ...
    tmp_a = phrase_trials(i,:,1);
    tmp_b = phrase_trials(i,:,2);
    tmp_c = phrase_trials(i,:,3);
    tmp_d = phrase_trials(i,:,4);
    
    %% PO units ...
    tmp_e = binding(tmp_a,dora_ps,1)'; tmp_e = eegfilt(tmp_e,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0);
    tmp_f = binding(tmp_b+tmp_c,dora_ps,2)'; tmp_f = eegfilt(tmp_f,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0);
    tmp_g = binding(tmp_d,dora_ps,1)'; tmp_g = eegfilt(tmp_g,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0);
    
    %% RB units ...
    tmp_h = binding(tmp_e,dora_ps,1)'; tmp_h = eegfilt(tmp_h,dora_ps.fs,0,dora_ps.order_filtering(2),0,[],0,'fir1',0);
    tmp_i = binding(tmp_f + tmp_g,dora_ps,2)'; tmp_i = eegfilt(tmp_i,dora_ps.fs,0,dora_ps.order_filtering(2),0,[],0,'fir1',0);
    
    %% P units ...
    tmp_j = binding(tmp_h + tmp_i,dora_ps,2)';tmp_j = eegfilt(tmp_j,dora_ps.fs,0,dora_ps.order_filtering(3),0,[],0,'fir1',0);
    
    %% organizing phrase results  ...
    phrase_units(strcmpi({phrase_units.id},'a')).response(i,:) = tmp_a./max(tmp_a);
    phrase_units(strcmpi({phrase_units.id},'b')).response(i,:) = tmp_b./max(tmp_b);
    phrase_units(strcmpi({phrase_units.id},'c')).response(i,:) = tmp_c./max(tmp_c);
    phrase_units(strcmpi({phrase_units.id},'d')).response(i,:) = tmp_d./max(tmp_d);
    phrase_units(strcmpi({phrase_units.id},'e')).response(i,:) = tmp_e./max(tmp_e);
    phrase_units(strcmpi({phrase_units.id},'f')).response(i,:) = tmp_f./max(tmp_f);
    phrase_units(strcmpi({phrase_units.id},'g')).response(i,:) = tmp_g./max(tmp_g);
    phrase_units(strcmpi({phrase_units.id},'h')).response(i,:) = tmp_h./max(tmp_h);
    phrase_units(strcmpi({phrase_units.id},'i')).response(i,:) = tmp_i./max(tmp_i);
    phrase_units(strcmpi({phrase_units.id},'j')).response(i,:) = tmp_j./max(tmp_j);
    
    fprintf([repmat('\b',1,n_str),'phrase simulations (%3d of %3d) ... '],i,dora_ps.n_trials);

end
fprintf('\n');
%% sentence simulation ...
n_str = fprintf('sentence simulations (%3d of %3d) ... ',0,0);

for i=1:dora_ps.n_trials
    %% syllable input ...
    tmp_a = sentence_trials(i,:,1);
    tmp_b = sentence_trials(i,:,2);
    tmp_c = sentence_trials(i,:,3);
    tmp_d = sentence_trials(i,:,4);
    
    %% PO units ...
    tmp_e = binding(tmp_a,dora_ps,1)'; tmp_e = eegfilt(tmp_e,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0);
    tmp_f = binding(tmp_b,dora_ps,1)'; tmp_f = eegfilt(tmp_f,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0);
    tmp_g = binding(tmp_c,dora_ps,3)'; tmp_g = eegfilt(tmp_g,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0); % 2 for beautiful plot ...
    tmp_h = binding(tmp_d,dora_ps,1)'; tmp_h = eegfilt(tmp_h,dora_ps.fs,0,dora_ps.order_filtering(1),0,[],0,'fir1',0);
    
    %% RB units ...
    tmp_i = binding(tmp_e + tmp_f,dora_ps,2)'; tmp_i = eegfilt(tmp_i,dora_ps.fs,0,dora_ps.order_filtering(2),0,[],0,'fir1',0);
    tmp_j = binding(tmp_g + tmp_h,dora_ps,2)'; tmp_j = eegfilt(tmp_j,dora_ps.fs,0,dora_ps.order_filtering(2),0,[],0,'fir1',0);
    
    %% P units ...
    tmp_k = binding(tmp_i + tmp_j,dora_ps,2)'; tmp_k = eegfilt(tmp_k,dora_ps.fs,0,dora_ps.order_filtering(3),0,[],0,'fir1',0);
    
    %% organizing sentence results  ...
    sentence_units(strcmpi({sentence_units.id},'a')).response(i,:) = tmp_a./max(tmp_a);
    sentence_units(strcmpi({sentence_units.id},'b')).response(i,:) = tmp_b./max(tmp_b);
    sentence_units(strcmpi({sentence_units.id},'c')).response(i,:) = tmp_c./max(tmp_c);
    sentence_units(strcmpi({sentence_units.id},'d')).response(i,:) = tmp_d./max(tmp_d);
    sentence_units(strcmpi({sentence_units.id},'e')).response(i,:) = tmp_e./max(tmp_e);
    sentence_units(strcmpi({sentence_units.id},'f')).response(i,:) = tmp_f./max(tmp_f);
    sentence_units(strcmpi({sentence_units.id},'g')).response(i,:) = tmp_g./max(tmp_g);
    sentence_units(strcmpi({sentence_units.id},'h')).response(i,:) = tmp_h./max(tmp_h);
    sentence_units(strcmpi({sentence_units.id},'i')).response(i,:) = tmp_i./max(tmp_i);
    sentence_units(strcmpi({sentence_units.id},'j')).response(i,:) = tmp_j./max(tmp_j);
    sentence_units(strcmpi({sentence_units.id},'k')).response(i,:) = tmp_k./max(tmp_k);
    
    fprintf([repmat('\b',1,n_str),'sentence simulations (%3d of %3d) ... '],i,dora_ps.n_trials);

end
fprintf('\n\n');
%% saving results ...
fprintf('saving simulation results ... \n');
save(fullfile(output_path,'simulation_results'),'phrase_units','sentence_units');

%% padding and binding ...
function data = binding(data,dora_ps,order_idx)
rng('shuffle');
n_smooth = dora_ps.n_binding(order_idx);
% sort_data = sort(data,'ascend');
median_sort_data = 0.90*median(data);
% baseline = repmat(median_sort_data,1,50);
pad_data = repmat(median_sort_data,1,n_smooth);
data = [pad_data,data,pad_data];
data = smooth(data,n_smooth);
data = data(n_smooth+1:end-n_smooth);

