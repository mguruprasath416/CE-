clc; clear; close all;

tau = 75e-6;
fs = 44100;
a = exp(-1/(tau*fs));
fprintf('Filter coefficient (a) = %.5f\n', a);

t = 0:1/fs:0.02;
x = sin(2*pi*1000*t) + 0.5*sin(2*pi*6000*t);
x = x / max(abs(x));

noise = 0.05*randn(size(x));
x_noisy = x + noise;

y_pre = filter([1 -a], 1, x_noisy);
y_de = filter(1, [1 -a], y_pre);

N = 2048;
[H_pre, f] = freqz([1 -a], 1, N, fs);
[H_de, ~] = freqz(1, [1 -a], N, fs);

figure;
subplot(2,1,1);
semilogx(f, 20*log10(abs(H_pre)));
title('Pre-Emphasis Filter Frequency Response');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
subplot(2,1,2);
semilogx(f, 20*log10(abs(H_de)));
title('De-Emphasis Filter Frequency Response');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;

figure;
subplot(3,1,1);
plot(t*1000, x_noisy); title('Original Noisy Signal');
xlabel('Time (ms)'); ylabel('Amplitude');
subplot(3,1,2);
plot(t*1000, y_pre); title('After Pre-Emphasis');
xlabel('Time (ms)'); ylabel('Amplitude');
subplot(3,1,3);
plot(t*1000, y_de); title('After De-Emphasis');
xlabel('Time (ms)'); ylabel('Amplitude');

SNR_in = 10*log10(var(x)/var(x_noisy - x));
SNR_out = 10*log10(var(x)/var(y_de - x));
fprintf('\nSNR before equalization = %.2f dB\n', SNR_in);
fprintf('SNR after equalization  = %.2f dB\n', SNR_out);
fprintf('SNR Improvement = %.2f dB\n', SNR_out - SNR_in);
