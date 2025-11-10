clc; clear; close all;

Am = 1;
Ac = 1;
fm = 1e3;
fc = 100e3;
mu = 0.7;
fs = 5e6;
t = 0:1/fs:2e-3;

m_t = Am * cos(2*pi*fm*t);
s_am = Ac * (1 + mu*m_t) .* cos(2*pi*fc*t);

SNR_levels = [10 20 30];
SNRout_vals = zeros(1,length(SNR_levels));

figure;
for i = 1:length(SNR_levels)
    SNRin_dB = SNR_levels(i);
    SNRin = 10^(SNRin_dB/10);
    noise = sqrt(var(s_am)/SNRin) * randn(size(s_am));
    s_noisy = s_am + noise;
    
    rectified = abs(s_noisy);
    N = 500;
    b = (1/N)*ones(1,N);
    m_rec = filter(b,1,rectified);
    m_rec = m_rec - mean(m_rec);
    
    m_aligned = m_t(1:length(m_rec));
    signal_power = var(m_aligned);
    noise_power_out = var(m_aligned - m_rec);
    SNRout = 10*log10(signal_power/noise_power_out);
    SNRout_vals(i) = SNRout;
    
    subplot(3,1,i);
    plot(t*1000, m_rec);
    title(['Recovered Signal at SNR = ', num2str(SNRin_dB), ' dB']);
    xlabel('Time (ms)'); ylabel('Amplitude');
end

figure;
plot(SNR_levels, SNRout_vals, '-o','LineWidth',2);
title('SNRout vs SNRin (AM System)');
xlabel('Input SNR (dB)');
ylabel('Output SNR (dB)');
grid on;

fprintf('\n--- AM System Noise Performance ---\n');
for i = 1:length(SNR_levels)
    fprintf('SNRin = %2d dB → SNRout = %.2f dB\n', SNR_levels(i), SNRout_vals(i));
end
