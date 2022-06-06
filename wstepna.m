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

%% Charakterystyka częstotliwościowa
k = 0:(N/2 - 1);
omega = 2 * pi * k / (N * Tp);

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
title('Charakterystyka amplitudowa');
xlabel('\omega [rad/s]');
grid on

%% Analiza korelacyjna - odpowiedź impulsowa
M = 30;
ryu = zeros(M, 1);
for tau=0:M-1
    ryu(tau + 1) = Covar([out, in], tau);
end

Ruu = zeros(M, M);
for i=0:M-1
    for j=0:M-1
        Ruu(i+1, j+1) = Covar([in, in], j - i);
    end
end

gM = 1 / Tp * pinv(Ruu) * ryu;
t = ((0:M-1) * Tp)';

fig_gM = figure('Position', [100 100 1200 900],...
                'Name', 'Odpowiedź impulsowa',...
                'NumberTitle', 'off');
plot(t, gM);
title('Odpowiedź impulsowa - analiza korelacyjna');
xlabel('t [s]');
grid on

%% Analiza korelacyjna - odpowiedź skokowa
hM = cumsum(gM * Tp);

fig_hM = figure('Position', [100 100 1200 900],...
                'Name', 'Odpowiedź skokowa',...
                'NumberTitle', 'off');
plot(t, hM);
title('Odpowiedź skokowa');
grid on

K = 0.96;
T3 = 0.48 / 2.67;
T4 = 0.48 / 3.67;

ag = 3;
tg = 0.38;
sg = 0.38;

sm = ag*t + sg - ag * tg;

s = tf('s');
Gm2 = K / (T3*s+1)^3;
Gm3 = K / (T4*s+1)^4;
T0 = (ag*tg - sg) / ag;
T = (K - sg + ag*tg) / ag - T0;
Gm0 = tf(K, [T, 1], 'IODelay', T0);

[y_2, ~] = step(Gm2, t);
[y_3, ~] = step(Gm3, t);
[y_0, ~] = step(Gm0, t);
hold on
plot(t, y_2);
plot(t, y_3);
plot(t, y_0);
legend('odp. analiza korelacyjna', 'model 3. rzędu', 'model 4. rzędu', 'inercja z opóźnieniem')
xlabel('t [s]')

% figure;
% plot(t, hM);
% hold on
% plot(t, sm)
% ylim([0, 13])
