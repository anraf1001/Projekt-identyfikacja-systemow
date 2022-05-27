close all; clear; clc;

%% Zainicjowanie zmiennych globalnych
Tp = 0.08;

%% Wczytanie oraz wyświetlenie danych pomiarowych
load('dane.mat');
N = length(in);
t = (0:length(in) - 1) * Tp;

fig_meas_data = figure('Position', [100 100 1200 900],...
                       'Name', 'Dane pomiarowe',...
                       'NumberTitle', 'off');

subplot(2, 1, 1);
plot(t, in);
title('Sygnał wejściowy - napięcie nagrzewnicy');
xlabel('t [s]');
grid on

subplot(2, 1, 2);
plot(t, out);
title('Sygnał wyjściowy - napięcie termopary');
xlabel('t [s]');
grid on

%% FFT sygnału wejściowego
k = 0:(N/2 - 1);
omega = 2 * pi * k / (N * Tp);

u_fft = fft(in);

fig_u_fft = figure('Position', [100 100 1200 900],...
                   'Name', 'FFT - sygnał wejściowy',...
                   'NumberTitle', 'off');

plot(omega, abs(u_fft(k + 1)));
title('FFT - sygnał wejściowy');
xlabel('\omega [rad/s]');
grid on

%% Charakterystyka częstotliwościowa
Mw = 150;
tau = -Mw:Mw;
wH = 0.5 * (1 + cos(tau * pi / Mw));

yu_corr = zeros(size(tau));
for i=1:length(tau)
    yu_corr(i) = Covar([out, in], tau(i));
end

uu_corr = zeros(size(tau));
for i=1:length(tau)
    uu_corr(i) = Covar([in, in], tau(i));
end

Phi_yu = zeros(1, N);
for i=0:N-1
    Phi_yu(i + 1) = Tp * sum(wH .* yu_corr .* exp(-1j .* tau .* (2*pi/N) * i));
end

Phi_uu = zeros(1, N);
for i=0:N-1
    Phi_uu(i + 1) = Tp * sum(wH .* uu_corr .* exp(-1j .* tau .* (2*pi/N) * i));
end

GN = Phi_yu ./ Phi_uu;

fig_bode = figure('Position', [100 100 1200 900],...
                  'Name', 'Charakterystyka częstotliwościowa',...
                  'NumberTitle', 'off');

semilogx(omega, 20 * log10(abs(GN(k + 1))));
title('Charakterystyka częstotliwościowa');
xlabel('\omega [rad/s]');
grid on
