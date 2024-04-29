function s03_envelope_info_extraction(paths,ps)
input_root_path = paths.stimuli;
output_root_path = paths.envelope_fft;

all_phrases = dir(fullfile(input_root_path,'phrases_exp','*.wav'));
all_sentences = dir(fullfile(input_root_path,'sentences_exp','*.wav'));

fs_downsample = ps.downsample;
%% downsample envelopes  ...
fft_envelope = [];
for audio_i = 1:length(all_phrases)
    %% reading files ...
    tmp_phrase_file_name = all_phrases(audio_i).name;
    tmp_sentence_file_name = all_sentences(audio_i).name;
    
    tmp_phrase = fullfile(input_root_path,'phrases_exp',tmp_phrase_file_name);
    tmp_sentence = fullfile(input_root_path,'sentences_exp',tmp_sentence_file_name);
    
    [x_phrase,~] = audioread(tmp_phrase);
    [x_sentence,fs] = audioread(tmp_sentence);
    
    %% envelope extraction ...
    phrase_envelope = abs(hilbert(x_phrase)); % sum(abs(y_phrase));
    %     phrase_envelope = eegfilt(phrase_envelope',fs,0,ps.lowpass_edge,0,[],[],'fir1');
    phrase_envelope = resample(phrase_envelope,ps.downsample,fs); % downsample ...
    
    sentence_envelope = abs(hilbert(x_sentence)); % sum(abs(y_sentence));
    %     sentence_envelope = eegfilt(sentence_envelope',fs,0,ps.lowpass_edge,0,[],[],'fir1');
    sentence_envelope = resample(sentence_envelope,ps.downsample,fs); % downsample ...
    
    
    %% plot fft envelope  ...
    signal_length = length(phrase_envelope);
    f = fs_downsample*(0:(signal_length/2))/signal_length;
    
    power_phrase = abs(fft(phrase_envelope)/length(phrase_envelope));
    power_phrase_single = power_phrase(1:length(phrase_envelope)/2+1);
    power_phrase_single(2:end-1) = 2*power_phrase_single(2:end-1);
    
    power_sentence = abs(fft(sentence_envelope)/length(sentence_envelope));
    power_sentence_single = power_sentence(1:length(sentence_envelope)/2+1);
    power_sentence_single(2:end-1) = 2*power_sentence_single(2:end-1);
    
    %% saving envelopes'frequency responses ...
    tmp_fft_envelope(1).name = tmp_phrase_file_name;
    tmp_fft_envelope(1).audio_type = 'phrase';
    tmp_fft_envelope(1).frequency = f(2:end)';
    tmp_fft_envelope(1).freq_response = power_phrase_single(2:end);
    tmp_fft_envelope(1).envelope = phrase_envelope;
    
    tmp_fft_envelope(2).name = tmp_sentence_file_name;
    tmp_fft_envelope(2).audio_type = 'sentence';
    tmp_fft_envelope(2).frequency = f(2:end)';
    tmp_fft_envelope(2).freq_response = power_sentence_single(2:end);
    tmp_fft_envelope(2).envelope = sentence_envelope;
    
    
    fft_envelope = [fft_envelope,tmp_fft_envelope];
end
save(fullfile(output_root_path,'envelope_fft.mat'),'fft_envelope');