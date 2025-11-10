clc; clear; close all;

fs = 8000;
fm = 1000;
t = 0:1/fs:5e-3;
Vm = 1;
x = Vm*sin(2*pi*fm*t);

n_bits = 4;
L = 2^n_bits;
xmax = Vm; xmin = -Vm;
delta_pcm = (xmax - xmin) / L;

xq = delta_pcm * floor(x/delta_pcm + 0.5);
q_noise_pcm = x - xq;
Pq_pcm = mean(q_noise_pcm.^2);
SNR_pcm = 10*log10(mean(x.^2)/Pq_pcm);

delta_dm = 0.1;
y_dm = zeros(size(x));
for i = 2:length(x)
    if x(i) >= y_dm(i-1)
        y_dm(i) = y_dm(i-1) + delta_dm;
    else
        y_dm(i) = y_dm(i-1) - delta_dm;
    end
end

q_noise_dm = x - y_dm;
Pq_dm = mean(q_noise_dm.^2);
SNR_dm = 10*log10(mean(x.^2)/Pq_dm);

figure;
subplot(3,1,1);
plot(t*1000, x, 'k', 'LineWidth', 1.5);
title('Original Message Signal');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
subplot(3,1,2);
plot(t*1000, xq, 'b', 'LineWidth', 1.2);
title('Reconstructed PCM Signal (4-bit Quantization)');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
subplot(3,1,3);
plot(t*1000, y_dm, 'r', 'LineWidth', 1.2);
title('Reconstructed Delta Modulated Signal');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;

fprintf('--- Quantization and SNR Analysis ---\n');
fprintf('PCM Quantization Step Size     = %.3f V\n', delta_pcm);
fprintf('PCM Quantization Noise Power   = %.6f\n', Pq_pcm);
fprintf('PCM SNR                        = %.2f dB\n', SNR_pcm);
fprintf('\nDelta Modulation Step Size     = %.3f V\n', delta_dm);
fprintf('DM Quantization Noise Power    = %.6f\n', Pq_dm);
fprintf('DM SNR                         = %.2f dB\n', SNR_dm);

figure;
plot(t*1000, x, 'k--', 'LineWidth', 1.2); hold on;
plot(t*1000, xq, 'b', 'LineWidth', 1.1);
plot(t*1000, y_dm, 'r', 'LineWidth', 1.1);
title('Comparison of PCM and DM Reconstructed Outputs');
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Original', 'PCM', 'Delta Modulation'); grid on;
