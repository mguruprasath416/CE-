clc; clear; close all;

fm = 1e3;
fs = 10e3;
fs_carrier = 100e3;
t = 0:1/fs_carrier:5e-3;

m_t = 0.5*(1 + sin(2*pi*fm*t));

n = 0:1/fs:5e-3;
m_samp = 0.5*(1 + sin(2*pi*fm*n));
pam_t = zeros(size(t));
for i = 1:length(n)
    [~, idx] = min(abs(t - n(i)));
    pam_t(idx) = m_samp(i);
end

triangular = sawtooth(2*pi*fs*t, 0.5);
pwm_t = double(m_t > (triangular + 1)/2);

ppd = 0.001;
ppm_t = zeros(size(t));
for i = 1:length(n)-1
    pulse_delay = m_samp(i) * (1/fs);
    start_time = n(i) + pulse_delay;
    end_time = start_time + (1/fs)/10;
    ppm_t(t >= start_time & t < end_time) = 1;
end

figure;
subplot(3,1,1);
plot(t*1000, pam_t, 'LineWidth', 1.2);
title('Pulse Amplitude Modulation (PAM)');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
subplot(3,1,2);
plot(t*1000, pwm_t, 'LineWidth', 1.2);
title('Pulse Width Modulation (PWM)');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
subplot(3,1,3);
plot(t*1000, ppm_t, 'LineWidth', 1.2);
title('Pulse Position Modulation (PPM)');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;

pam_amp = max(pam_t);
avg_pwm_width = mean(diff(find(diff(pwm_t)==1)))/fs_carrier;
avg_ppm_shift = mean(m_samp)/fs;

fprintf('\n--- PULSE PARAMETERS ---\n');
fprintf('PAM max amplitude        = %.2f\n', pam_amp);
fprintf('Average PWM pulse width  = %.6f s\n', avg_pwm_width);
fprintf('Average PPM time shift   = %.6f s\n', avg_ppm_shift);

N = length(t);
f = (-N/2:N/2-1)*(fs_carrier/N);
PAM_spec = abs(fftshift(fft(pam_t)))/N;
PWM_spec = abs(fftshift(fft(pwm_t)))/N;
PPM_spec = abs(fftshift(fft(ppm_t)))/N;

figure;
subplot(3,1,1); plot(f/1000, PAM_spec); title('PAM Spectrum'); xlim([0 50]);
xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); grid on;
subplot(3,1,2); plot(f/1000, PWM_spec); title('PWM Spectrum'); xlim([0 50]);
xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); grid on;
subplot(3,1,3); plot(f/1000, PPM_spec); title('PPM Spectrum'); xlim([0 50]);
xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); grid on;

noise_tolerance_PAM = 1 / var(pam_t(pam_t>0));
noise_tolerance_PWM = 1 / var(pwm_t(pwm_t>0));
noise_tolerance_PPM = 1 / var(ppm_t(ppm_t>0));

fprintf('\n--- NOISE TOLERANCE (higher = better) ---\n');
fprintf('PAM  : %.3f\n', noise_tolerance_PAM);
fprintf('PWM  : %.3f\n', noise_tolerance_PWM);
fprintf('PPM  : %.3f\n', noise_tolerance_PPM);

bw_pam = 10e3;
bw_pwm = 15e3;
bw_ppm = 12e3;
bit_eff_pam = 1/bw_pam;
bit_eff_pwm = 1/bw_pwm;
bit_eff_ppm = 1/bw_ppm;

fprintf('\n--- SPECTRAL EFFICIENCY (1/BW) ---\n');
fprintf('PAM  : %.3e\n', bit_eff_pam);
fprintf('PWM  : %.3e\n', bit_eff_pwm);
fprintf('PPM  : %.3e\n', bit_eff_ppm);
