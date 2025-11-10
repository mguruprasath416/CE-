clc; clear; close all;

fs = 5e6;
t = 0:1/fs:2e-3;
fm = 5e3;
fc = 1e6;
fIF = 455e3;
fLO = fc - fIF;
Ac = 1; Am = 0.5;
mu = Am/Ac;

m_t = Am*cos(2*pi*fm*t);
s_RF = Ac*(1 + mu*cos(2*pi*fm*t)).*cos(2*pi*fc*t);

LO = cos(2*pi*fLO*t);
s_mixed = s_RF .* LO;

bpFilt = fir1(500, [(fIF-10e3)/(fs/2) (fIF+10e3)/(fs/2)], 'bandpass');
s_IF = filter(bpFilt, 1, s_mixed);

s_rectified = abs(s_IF);
lpf = fir1(200, (20e3)/(fs/2));
s_baseband = filter(lpf, 1, s_rectified);
s_baseband = s_baseband - mean(s_baseband);

N = length(t);
f = (-N/2:N/2-1)*(fs/N);
RF_spec = abs(fftshift(fft(s_RF)))/N;
IF_spec = abs(fftshift(fft(s_IF)))/N;
BB_spec = abs(fftshift(fft(s_baseband)))/N;

figure;
subplot(3,1,1);
plot(f/1e3, RF_spec); xlim([0 2000]);
title('RF Spectrum (AM at 1 MHz)'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|');
subplot(3,1,2);
plot(f/1e3, IF_spec); xlim([0 1000]);
title('IF Spectrum (Centered at 455 kHz)'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|');
subplot(3,1,3);
plot(f/1e3, BB_spec); xlim([0 50]);
title('Baseband Message Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|');

figure;
subplot(3,1,1); plot(t(1:2000)*1e3, s_RF(1:2000)); title('RF AM Signal'); xlabel('Time (ms)');
subplot(3,1,2); plot(t(1:2000)*1e3, s_IF(1:2000)); title('IF Signal (455 kHz)'); xlabel('Time (ms)');
subplot(3,1,3); plot(t(1:2000)*1e3, s_baseband(1:2000)); title('Recovered Baseband Message'); xlabel('Time (ms)');

fprintf('--- SUPERHETERODYNE RECEIVER RESULTS ---\n');
fprintf('RF Carrier Frequency   = %.0f kHz\n', fc/1e3);
fprintf('LO Frequency           = %.0f kHz\n', fLO/1e3);
fprintf('Intermediate Frequency = %.0f kHz\n', fIF/1e3);
fprintf('Baseband Message Freq  = %.0f kHz\n', fm/1e3);
