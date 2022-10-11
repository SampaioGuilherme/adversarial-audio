clear;

sinal = 'forecast.wav';

[sig1, fs1] = audioread(sinal); % Importa arquivo de audio
sig1 = sig1(:, 1); % Audio em mono, conforme gravacao

fs = 20000; % Define taxa de amostragem a 20 kHz
x = resample(sig1, fs, fs1); % Faz resample do sinal original

frm_len = 0.02*fs; % Define tamanho da janela de 20 ms
frm_sft = 0.01*fs; % Deslocamento de 10 ms, sobreposicao de 50%

remainder = mod(length(x), frm_sft);
x_trim = x(1:(end-remainder)); % Normaliza o tamanho do sinal
n_frames = length(x_trim)/frm_sft; % Obtem o numero de janelas

X_sft = reshape(x_trim, [frm_sft, n_frames]); % Janela o sinal
X_frames = [zeros(frm_sft,1) X_sft(:,1:(end-1)); X_sft]; % Monta janelas
X = X_frames.*hamming(frm_len); 

X_fft = fft(X); % Realiza a FFT do sinal
X_mag = abs(X_fft);
X_phase = angle(X_fft);

% X_fft_rec = ifft(X_fft);
X_fft_rec = X_mag.*exp(1i*X_phase); % Reconstroi o sinal
x_rec = real(ifft(X_fft_rec)./hamming(frm_len));
x_rec_aux = x_rec((frm_len/4+1):(frm_len*3/4),:);
% x_rec_aux = X_rec(1:(frm_len/2),:);
% x_rec_aux = X_rec((frm_len/2+1):end ,:);

sig_rec = reshape(x_rec_aux, [length(x_trim), 1]);

subplot(2,1,1)
plot(real(x_rec_aux(:)))
title('Sinal reconstruído')

rand_phase = rand(frm_len/2+1,n_frames)*2*pi - pi;
X_rand_phase = [rand_phase; -flipud(rand_phase(2:(frm_len/2),:))];
x_rec_rand = real(ifft(X_mag.*exp(1i*X_rand_phase))./hamming(frm_len));
x_rec_aux_rand = x_rec_rand((frm_len/4+1):(frm_len*3/4),:); 

subplot(2,1,2)
plot(x_rec_aux_rand(:))
title('Sinal reconstruído a partir de fase aleatória')

audiowrite(append('rec_', sinal), x_rec_aux(:), fs);
audiowrite(append('rand_phase_rec_', sinal), real(x_rec_aux_rand(:)), fs);

