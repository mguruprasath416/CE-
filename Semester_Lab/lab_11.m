clc; clear; close all;

Am = 1;
fm = 1000;
fs = 8000;
t = 0:1/fs:3e-3;
L = 8;
n = log2(L);

x = Am * sin(2*pi*fm*t);

xmin = -Am; xmax = Am;
delta = (xmax - xmin)/L;
partition = xmin+delta:delta:xmax-delta;
codebook = xmin+delta/2:delta:xmax-delta/2;

index = zeros(size(x));
quantized = zeros(size(x));
for i = 1:length(x)
    idx = find(x(i) < partition, 1);
    if isempty(idx), idx = L; end
    index(i) = idx;
    quantized(i) = codebook(idx);
end

bit_stream = [];
for i = 1:length(index)
    bits = dec2bin(index(i)-1, n) - '0';
    bit_stream = [bit_stream bits];
end

figure;
subplot(3,1,1); plot(t*1000, x, 'b');
title('Original Message Signal'); xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
subplot(3,1,2); stem(t*1000, quantized, 'filled');
title('Quantized Signal (8 Levels)'); xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
subplot(3,1,3); stairs(bit_stream(1:60), 'LineWidth', 1.5);
title('PCM Bitstream (first 60 bits)'); ylabel('Bit'); xlabel('Sample Index'); grid on;

N = length(bit_stream);
f = (-N/2:N/2-1)*(fs*N/length(x))/N;
bit_fft = abs(fftshift(fft(bit_stream)))/N;

figure;
plot(f/1000, bit_fft);
title('Spectrum of PCM Bitstream'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|');
xlim([0 50]); grid on;

PSD = (bit_fft.^2)/N;
figure;
plot(f/1000, PSD);
title('Power Spectral Density of PCM Signal'); xlabel('Frequency (kHz)'); ylabel('Power');
xlim([0 50]); grid on;

lpf_cutoff = 1500;
b = fir1(128, lpf_cutoff/(fs/2));
x_rec = filter(b, 1, quantized);

figure;
plot(t*1000, x, 'b', 'LineWidth', 1.2); hold on;
plot(t*1000, x_rec, 'r--', 'LineWidth', 1.2);
title('Original vs Reconstructed Signal');
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Original','Reconstructed'); grid on;

MSE = mean((x - x_rec).^2);
BW = fs * n / 2;

fprintf('--- PCM SYSTEM RESULTS ---\n');
fprintf('Quantization Levels (L) = %d\n', L);
fprintf('Bits per Sample (n)     = %.0f\n', n);
fprintf('Quantization Step Δ     = %.4f\n', delta);
fprintf('Reconstruction MSE      = %.6f\n', MSE);
fprintf('Estimated PCM Bandwidth = %.1f Hz\n', BW);
