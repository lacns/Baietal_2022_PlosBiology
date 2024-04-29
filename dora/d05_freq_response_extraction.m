function d05_freq_response_extraction(paths,dora_ps)
input_path = paths.result_path;
output_path = input_path;
load(fullfile(input_path,'simulation_results'));


%% frequencies ...
sig_length = dora_ps.fs*dora_ps.trial_length;
f = dora_ps.fs.*(0:(sig_length/2))/sig_length;
f=f(2:end);


%% phrases ...
n_str = fprintf('phrase frequency decomposition (%3d of %3d)...',0,0);
for i=1:length(phrase_units)
    tmp_unit_id = phrase_units(i).id;
    tmp_unit_response = phrase_units(i).response';
    
    %% fft unit activation ...
    tmp_unit_fft = fft(tmp_unit_response); %
    tmp_unit_fft = tmp_unit_fft(2:length(f)+1,:);
    
    %% power extraction ...
    tmp_unit_power = 2*(abs(tmp_unit_fft/sig_length));
    
    %% phase coherence ...
    tmp_unit_itpc = phase_coherence(tmp_unit_fft,dora_ps);
    
    %% organizing results ...
    phrase_units(i).frex = f;
    phrase_units(i).power = tmp_unit_power;
    phrase_units(i).itpc = tmp_unit_itpc;
    phrase_units(i).complex = tmp_unit_fft;
    
    fprintf([repmat('\b',1,n_str),'phrase frequency decomposition (%3d of %3d)...'],i,length(phrase_units));
    
end

fprintf('\n');

%% sentences ...
n_str = fprintf('sentence frequency decomposition (%3d of %3d)...',0,0);

for i=1:length(sentence_units)
    tmp_unit_id = sentence_units(i).id;
    tmp_unit_response = sentence_units(i).response';
    
    %% fft unit activation ...
    tmp_unit_fft = fft(tmp_unit_response); %
    tmp_unit_fft = tmp_unit_fft(2:length(f)+1,:);
    
    %% power extraction ...
    tmp_unit_power = 2*(abs(tmp_unit_fft/sig_length));
    
    %% phase coherence ...
    tmp_unit_itpc = phase_coherence(tmp_unit_fft,dora_ps);
    
    %% organizing results ...
    sentence_units(i).frex = f;
    sentence_units(i).power = tmp_unit_power;
    sentence_units(i).itpc = tmp_unit_itpc;
    sentence_units(i).complex = tmp_unit_fft;
    
    fprintf([repmat('\b',1,n_str),'sentence frequency decomposition (%3d of %3d)...'],i,length(sentence_units));
    
end
fprintf('\n');
fprintf('saving frequency decomposition results ...\n\n');
save(fullfile(output_path,'power_phase'),'phrase_units','sentence_units');

function itpc = phase_coherence(data,dora_ps)
rng('shuffle');
for i=1:dora_ps.n_itpc
    tmp_data = data(:,randperm(size(data,2),dora_ps.n_per_calculation));
    itpc(:,i) = abs(mean(exp(1i.*angle(tmp_data)),2));
end
