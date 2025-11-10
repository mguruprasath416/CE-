clc; clear; close all;

N = 1e4;
Rb = 1000;
fs = 10*Rb;
Tb = 1/Rb;
t = 0:1/fs:Tb-1/fs;
A = 1;
fc = 2e3;
f1 = 2e3; f2 = 4e3;

bits = randi([0 1], 1, N);

ASK_sig = []; FSK_sig = []; BPSK_sig = [];
for i = 1:N
    if bits(i) == 1
        ASK_sig = [ASK_sig A*cos(2*pi*fc*t)];
        FSK_sig = [FSK_sig A*cos(2*pi*f1*t)];
        BPSK_sig = [BPSK_sig A*cos(2*pi*fc*t)];
    else
        ASK_sig = [ASK_sig zeros(1,length(t))];
        FSK_sig = [FSK_sig A*cos(2*pi*f2*t)];
        BPSK_sig = [BPSK_sig -A*cos(2*pi*fc*t)];
    end
end

SNR_dB = 0:2:15;
BER_ASK = zeros(size(SNR_dB));
BER_FSK = zeros(size(SNR_dB));
BER_BPSK = zeros(size(SNR_dB));

for k = 1:length(SNR_dB)
    SNR_lin = 10^(SNR_dB(k)/10);
    noise_std = sqrt(A^2/(2*SNR_lin));
    ASK_noisy = ASK_sig + noise_std*randn(size(ASK_sig));
    FSK_noisy = FSK_sig + noise_std*randn(size(FSK_sig));
    BPSK_noisy = BPSK_sig + noise_std*randn(size(BPSK_sig));

    rec_bits_ASK = zeros(1,N);
    rec_bits_FSK = zeros(1,N);
    rec_bits_BPSK = zeros(1,N);

    for i = 1:N
        idx = (i-1)*length(t) + (1:length(t));
        carrier = cos(2*pi*fc*t);
        rec_ask = trapz(t, ASK_noisy(idx).*carrier);
        rec_bits_ASK(i) = rec_ask > 0.3;
        rec_f1 = trapz(t, FSK_noisy(idx).*cos(2*pi*f1*t));
        rec_f2 = trapz(t, FSK_noisy(idx).*cos(2*pi*f2*t));
        rec_bits_FSK(i) = rec_f1 > rec_f2;
        rec_bpsk = trapz(t, BPSK_noisy(idx).*carrier);
        rec_bits_BPSK(i) = rec_bpsk > 0;
    end

    BER_ASK(k) = sum(bits ~= rec_bits_ASK)/N;
    BER_FSK(k) = sum(bits ~= rec_bits_FSK)/N;
    BER_BPSK(k) = sum(bits ~= rec_bits_BPSK)/N;
    fprintf('SNR = %2d dB: BER(ASK)=%.4f, BER(FSK)=%.4f, BER(BPSK)=%.4f\n', SNR_dB(k), BER_ASK(k), BER_FSK(k), BER_BPSK(k));
end

figure;
semilogy(SNR_dB, BER_ASK, '-or','LineWidth',1.5); hold on;
semilogy(SNR_dB, BER_FSK, '-sb','LineWidth',1.5);
semilogy(SNR_dB, BER_BPSK, '-^g','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title('BER Performance of ASK, FSK, and BPSK over AWGN');
legend('ASK','FSK','BPSK','Location','southwest');

figure;
subplot(3,1,1); plot(ASK_sig(1:10*length(t)),'r'); title('ASK Modulated Signal'); ylabel('Amplitude');
subplot(3,1,2); plot(FSK_sig(1:10*length(t)),'b'); title('FSK Modulated Signal'); ylabel('Amplitude');
subplot(3,1,3); plot(BPSK_sig(1:10*length(t)),'g'); title('BPSK Modulated Signal'); ylabel('Amplitude');
xlabel('Sample Index');
