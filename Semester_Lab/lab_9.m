clc; clear; close all;

fs = 1e5;
t = 0:1/fs:2e-3;
fm = 1e3;
fc = 10e3;
Ac = 1;
Am = 1;

m_t = Am*cos(2*pi*fm*t);

mu = 0.5;
s_am = Ac*(1 + mu*m_t).*cos(2*pi*fc*t);

df = 5e3;
beta = df/fm;
s_fm = Ac*cos(2*pi*fc*t + beta*sin(2*pi*fm*t));

Rb = 1e3;
bits = randi([0 1], 1, length(t));
s_ask = Ac*(bits).*cos(2*pi*fc*t);

bipolar = 2*bits - 1;
s_bpsk = Ac*bipolar.*cos(2*pi*fc*t);

SNR_dB = 15;
SNR_lin = 10^(SNR_dB/10);
noise_am = sqrt(var(s_am)/SNR_lin)*randn(size(s_am));
noise_fm = sqrt(var(s_fm)/SNR_lin)*randn(size(s_fm));
noise_ask = sqrt(var(s_ask)/SNR_lin)*randn(size(s_ask));
noise_bpsk = sqrt(var(s_bpsk)/SNR_lin)*randn(size(s_bpsk));

s_am_noisy = s_am + noise_am;
s_fm_noisy = s_fm + noise_fm;
s_ask_noisy = s_ask + noise_ask;
s_bpsk_noisy = s_bpsk + noise_bpsk;

N = length(t);
f = (-N/2:N/2-1)*(fs/N);
fft_am = abs(fftshift(fft(s_am_noisy)))/N;
fft_fm = abs(fftshift(fft(s_fm_noisy)))/N;
fft_ask = abs(fftshift(fft(s_ask_noisy)))/N;
fft_bpsk = abs(fftshift(fft(s_bpsk_noisy)))/N;

bw = @(s_fft,freqs) ...
     (freqs(find(cumsum(s_fft.^2)/sum(s_fft.^2)>0.005,1)) - ...
      freqs(find(cumsum(s_fft.^2)/sum(s_fft.^2)>0.995,1)));

BW_AM = abs(bw(fft_am,f));
BW_FM = abs(bw(fft_fm,f));
BW_ASK = abs(bw(fft_ask,f));
BW_BPSK = abs(bw(fft_bpsk,f));

calcSNR = @(sig,noisy) 10*log10(var(sig)/var(noisy-sig));

SNRout_AM = calcSNR(s_am,s_am_noisy);
SNRout_FM = calcSNR(s_fm,s_fm_noisy);
SNRout_ASK = calcSNR(s_ask,s_ask_noisy);
SNRout_BPSK = calcSNR(s_bpsk,s_bpsk_noisy);

eta_AM = (mu^2)/(2+mu^2);
eta_FM = 1/(1+beta^2);
eta_ASK = mean(bits.^2)/(mean(bits.^2)+var(noise_ask));
eta_BPSK = mean(bipolar.^2)/(mean(bipolar.^2)+var(noise_bpsk));

fprintf('\n--- COMPARATIVE RESULTS ---\n');
fprintf('Scheme\tBW(Hz)\tSNRout(dB)\tEfficiency\n');
fprintf('AM\t%.0f\t%.2f\t\t%.2f\n', BW_AM, SNRout_AM, eta_AM);
fprintf('FM\t%.0f\t%.2f\t\t%.2f\n', BW_FM, SNRout_FM, eta_FM);
fprintf('ASK\t%.0f\t%.2f\t\t%.2f\n', BW_ASK, SNRout_ASK, eta_ASK);
fprintf('BPSK\t%.0f\t%.2f\t\t%.2f\n', BW_BPSK, SNRout_BPSK, eta_BPSK);

schemes = {'AM','FM','ASK','BPSK'};
SNR_vals = [SNRout_AM SNRout_FM SNRout_ASK SNRout_BPSK];
Eff_vals = [eta_AM eta_FM eta_ASK eta_BPSK];

figure;
scatter(Eff_vals*100, SNR_vals, 100, 'filled');
text(Eff_vals*100+0.5, SNR_vals, schemes);
grid on;
xlabel('Power Efficiency (%)');
ylabel('Output SNR (dB)');
title('SNR vs Power Efficiency Comparison');

figure;
subplot(4,1,1); plot(t*1000, s_am_noisy); title('AM Noisy Signal'); ylabel('Amplitude');
subplot(4,1,2); plot(t*1000, s_fm_noisy); title('FM Noisy Signal'); ylabel('Amplitude');
subplot(4,1,3); plot(t*1000, s_ask_noisy); title('ASK Noisy Signal'); ylabel('Amplitude');
subplot(4,1,4); plot(t*1000, s_bpsk_noisy); title('BPSK Noisy Signal'); ylabel('Amplitude');
xlabel('Time (ms)');
