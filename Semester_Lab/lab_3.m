clc; clear; close all;

Ac = 1;
Am = 1;
fm = 1e3;
fc = 50e3;
fs = 1e6;
t = 0:1/fs:3e-3;

df_values = [1e3 3e3 5e3];

figure;
for i = 1:length(df_values)
    df = df_values(i);
    beta = df / fm;
    s_fm = Ac * cos(2*pi*fc*t + beta * sin(2*pi*fm*t));

    N = length(t);
    f = (-N/2:N/2-1)*(fs/N);
    S_FM = abs(fftshift(fft(s_fm)))/N;

    subplot(length(df_values),1,i);
    plot(f/1000, S_FM);
    title(['FM Spectrum (Δf = ', num2str(df/1e3), ' kHz, β = ', num2str(beta), ')']);
    xlabel('Frequency (kHz)'); ylabel('|Amplitude|');
    xlim([30 70]); grid on;

    BW_theoretical = 2*(df + fm);
    fprintf('Δf = %.0f Hz → β = %.2f → Carson Bandwidth = %.0f Hz\n', df, beta, BW_theoretical);
end
