clc; clear; close all;

Am = 1;
Ac = 1;
fm = 1e3;
fc = 10e3;
fs = 1e6;
t = 0:1/fs:2e-3;

m_t = Am * cos(2*pi*fm*t);
c_t = Ac * cos(2*pi*fc*t);

mu = 0.5;
s_am = Ac * (1 + mu * cos(2*pi*fm*t)) .* cos(2*pi*fc*t);

df = 5e3;
beta = df / fm;
s_fm = Ac * cos(2*pi*fc*t + beta * sin(2*pi*fm*t));

N = length(t);
f = (-N/2:N/2-1)*(fs/N);

AM_FFT = abs(fftshift(fft(s_am)))/N;
FM_FFT = abs(fftshift(fft(s_fm)))/N;

figure;
subplot(2,1,1); plot(t*1000, s_am); title('AM Signal (Time Domain)');
xlabel('Time (ms)'); ylabel('Amplitude');
subplot(2,1,2); plot(t*1000, s_fm); title('FM Signal (Time Domain)');
xlabel('Time (ms)'); ylabel('Amplitude');

figure;
subplot(2,1,1); plot(f/1000, AM_FFT);
title('AM Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|');
xlim([0 30]);
subplot(2,1,2); plot(f/1000, FM_FFT);
title('FM Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|');
xlim([0 30]);

AM_eff = (mu^2) / (2 + mu^2);
AM_BW = 2 * fm;
FM_BW = 2 * (df + fm);

fprintf('\n--- PARAMETERS ---\n');
fprintf('AM Modulation Index (mu) = %.2f\n', mu);
fprintf('AM Power Efficiency = %.2f%%\n', AM_eff*100);
fprintf('AM Bandwidth = %.1f Hz\n', AM_BW);
fprintf('\nFM Modulation Index (beta) = %.2f\n', beta);
fprintf('FM Bandwidth (Carson rule) = %.1f Hz\n', FM_BW);

SNR_dB = 20;
SNR = 10^(SNR_dB/10);
noise_am = sqrt(var(s_am)/SNR)*randn(size(s_am));
noise_fm = sqrt(var(s_fm)/SNR)*randn(size(s_fm));
s_am_noisy = s_am + noise_am;
s_fm_noisy = s_fm + noise_fm;

figure;
subplot(2,1,1);
plot(t*1000, s_am_noisy);
title('AM with Noise'); xlabel('Time (ms)'); ylabel('Amplitude');
subplot(2,1,2);
plot(t*1000, s_fm_noisy);
title('FM with Noise'); xlabel('Time (ms)'); ylabel('Amplitude');

SNR_AM_measured = 10*log10(var(s_am)/var(noise_am));
SNR_FM_measured = 10*log10(var(s_fm)/var(noise_fm));

fprintf('\nMeasured SNR (AM) = %.2f dB\n', SNR_AM_measured);
fprintf('Measured SNR (FM) = %.2f dB\n', SNR_FM_measured);
