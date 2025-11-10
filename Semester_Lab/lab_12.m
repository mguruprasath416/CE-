clc; clear; close all;

fs = 1e5;
t = 0:1/fs:3e-3;
fm = 1e3;
fc = 10e3;
Ac = 1; Am = 1;
SNRin_dB = 10;
SNRlin = 10^(SNRin_dB/10);
m_t = Am*cos(2*pi*fm*t);

mu = 0.7;
s_AM = Ac*(1 + mu*m_t).*cos(2*pi*fc*t);
noise_am = sqrt(var(s_AM)/SNRlin)*randn(size(s_AM));
r_AM = s_AM + noise_am;
env = abs(r_AM);
b = fir1(100, (2*fm)/(fs/2));
rec_AM = filter(b,1,env) - mean(env);
rec_AM = rec_AM(1:length(m_t));
SNRout_AM = 10*log10(var(m_t)/var(m_t - rec_AM));

kf = 2*pi*5e3;
s_FM = Ac*cos(2*pi*fc*t + kf*cumsum(m_t)/fs);
noise_fm = sqrt(var(s_FM)/SNRlin)*randn(size(s_FM));
r_FM = s_FM + noise_fm;
analytic_FM = hilbert(r_FM);
demod_FM = diff(unwrap(angle(analytic_FM)));
demod_FM = demod_FM - mean(demod_FM);
rec_FM = filter(b,1,demod_FM);
L = min(length(m_t), length(rec_FM));
rec_FM = rec_FM(1:L);
m_t_FM = m_t(1:L);
SNRout_FM = 10*log10(var(m_t_FM)/var(m_t_FM - rec_FM));

bits = randi([0 1], 1, length(t));
s_ASK = Ac*(bits).*cos(2*pi*fc*t);
noise_ask = sqrt(var(s_ASK)/SNRlin)*randn(size(s_ASK));
r_ASK = s_ASK + noise_ask;
rec_ask = r_ASK .* cos(2*pi*fc*t);
rec_bits = rec_ask > 0.3;
BER_ASK = sum(bits ~= rec_bits)/length(bits);
SNRout_ASK = 10*log10(var(bits)/var(bits - rec_bits));

N = 1000;
bits_qpsk = randi([0 1], 1, N);
M = 4; k = log2(M);
numSymbols = N/k;
symbols = zeros(1,numSymbols);
for i = 1:numSymbols
    b1 = bits_qpsk(2*i-1); b2 = bits_qpsk(2*i);
    if     b1==0 && b2==0, symbols(i) = (1 + 1j);
    elseif b1==0 && b2==1, symbols(i) = (-1 + 1j);
    elseif b1==1 && b2==1, symbols(i) = (-1 - 1j);
    else                   symbols(i) = (1 - 1j);
    end
end
symbols = symbols / sqrt(2);
samplesPerSymbol = 20;
s_QPSK = real(repelem(symbols, samplesPerSymbol));
noise_qpsk = sqrt(var(s_QPSK)/SNRlin)*randn(size(s_QPSK));
r_QPSK = s_QPSK + noise_qpsk;
rec_QPSK = sign(r_QPSK);
BER_QPSK = sum(sign(s_QPSK) ~= sign(r_QPSK)) / length(s_QPSK);
SNRout_QPSK = 10*log10(var(s_QPSK)/var(s_QPSK - r_QPSK));

Nfft = length(t);
f = (-Nfft/2:Nfft/2-1)*(fs/Nfft);
figure;
subplot(4,1,1); plot(f/1000, abs(fftshift(fft(s_AM)))/Nfft);
title('AM Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); xlim([0 30]);
subplot(4,1,2); plot(f/1000, abs(fftshift(fft(s_FM)))/Nfft);
title('FM Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); xlim([0 30]);
subplot(4,1,3); plot(f/1000, abs(fftshift(fft(s_ASK)))/Nfft);
title('ASK Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); xlim([0 30]);
subplot(4,1,4); plot(f/1000, abs(fftshift(fft(s_QPSK(1:Nfft))))/Nfft);
title('QPSK Spectrum'); xlabel('Frequency (kHz)'); ylabel('|Amplitude|'); xlim([0 30]);
sgtitle('Frequency Spectra of All Modulation Schemes');

schemes = {'AM','FM','ASK','QPSK'};
SNRout_vals = [SNRout_AM SNRout_FM SNRout_ASK SNRout_QPSK];
BER_vals = [NaN NaN BER_ASK BER_QPSK];
figure;
yyaxis left
bar(SNRout_vals); ylabel('SNR_{out} (dB)');
yyaxis right
plot(3:4, BER_vals(3:4), 'ro--', 'LineWidth', 1.5);
xticks(1:4); xticklabels(schemes);
title('SNR_{out} and BER Comparison Across Modulation Schemes');
xlabel('Modulation Type'); grid on;
legend('SNR_{out}','BER');
