clc; clear; close all;

N = 1000;
bits = randi([0 1], 1, N);
EbN0_dB = 10;
M = 4;
k = log2(M);
numSymbols = N/k;

symbols = zeros(1, numSymbols);
for i = 1:numSymbols
    b1 = bits(2*i-1);
    b2 = bits(2*i);
    if     b1==0 && b2==0, symbols(i) =  (1 + 1j);
    elseif b1==0 && b2==1, symbols(i) = (-1 + 1j);
    elseif b1==1 && b2==1, symbols(i) = (-1 - 1j);
    else                   symbols(i) = (1 - 1j);
    end
end
symbols = symbols / sqrt(2);

Es = mean(abs(symbols).^2);
EbN0 = 10^(EbN0_dB/10);
N0 = Es/(2*EbN0);
noise = sqrt(N0)*(randn(size(symbols)) + 1j*randn(size(symbols)));
rx = symbols + noise;

rx_bits = zeros(1, N);
for i = 1:numSymbols
    rI = real(rx(i));
    rQ = imag(rx(i));
    if     (rI>0 && rQ>0), rx_bits(2*i-1:2*i) = [0 0];
    elseif (rI<0 && rQ>0), rx_bits(2*i-1:2*i) = [0 1];
    elseif (rI<0 && rQ<0), rx_bits(2*i-1:2*i) = [1 1];
    else                   rx_bits(2*i-1:2*i) = [1 0];
    end
end

BER = sum(bits ~= rx_bits)/N;
fprintf('SNR = %d dB --> BER = %.5f\n', EbN0_dB, BER);

figure;
subplot(1,2,1);
plot(real(symbols), imag(symbols), 'bo'); axis equal; grid on;
title('Constellation Before Noise');
xlabel('In-phase'); ylabel('Quadrature');
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
subplot(1,2,2);
plot(real(rx), imag(rx), 'r.'); axis equal; grid on;
title(['Constellation After Noise (SNR = ', num2str(EbN0_dB), ' dB)']);
xlabel('In-phase'); ylabel('Quadrature');
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

samplesPerSymbol = 20;
t = 0:1/samplesPerSymbol:numSymbols-1/samplesPerSymbol;
tx_wave = repelem(real(symbols), samplesPerSymbol);
rx_wave = repelem(real(rx), samplesPerSymbol);
figure;
plot(reshape(rx_wave(1:400), samplesPerSymbol, []));
title('Eye Diagram (I-Channel)');
xlabel('Sample'); ylabel('Amplitude'); grid on;

SNR_dB = 0:2:15;
BER_sim = zeros(size(SNR_dB));
for k = 1:length(SNR_dB)
    EbN0 = 10^(SNR_dB(k)/10);
    N0 = Es/(2*EbN0);
    noise = sqrt(N0)*(randn(size(symbols)) + 1j*randn(size(symbols)));
    rx_temp = symbols + noise;
    rx_bits_temp = zeros(1, N);
    for i = 1:numSymbols
        rI = real(rx_temp(i)); rQ = imag(rx_temp(i));
        if     (rI>0 && rQ>0), rx_bits_temp(2*i-1:2*i) = [0 0];
        elseif (rI<0 && rQ>0), rx_bits_temp(2*i-1:2*i) = [0 1];
        elseif (rI<0 && rQ<0), rx_bits_temp(2*i-1:2*i) = [1 1];
        else                   rx_bits_temp(2*i-1:2*i) = [1 0];
        end
    end
    BER_sim(k) = sum(bits ~= rx_bits_temp)/N;
end

figure;
semilogy(SNR_dB, BER_sim, '-ob', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('QPSK BER vs SNR (Simulated)');
