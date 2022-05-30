close all; clear; clc;

%% Zainicjowanie zmiennych globalnych
Tp = 0.08;

%% Wczytanie danych pomiarowych
load('dane.mat');
N = length(in);
t = (0:length(in) - 1) * Tp;

%% Podział danych pomiarowych na 2 podzbiory
div = fix(N / 4);

in_est = in(1:div);
in_wer = in(div+1:end);

out_est = out(1:div);
out_wer = out(div+1:end);

N_est = length(in_est);
N_wer = length(in_wer);

%% LS - 3 rząd
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) u(n-1) u(n-2) u(n-3)]
Phi3 = zeros(N_est, 6);
Phi3(2, :) = [-out_est(1), 0, 0, in_est(1), 0, 0];
Phi3(3, :) = [-out_est(2), -out_est(1), 0, in_est(2), in_est(1), 0];

for k=4:N_est
    Phi3(k, :) = [-out_est(k-1), -out_est(k-2), -out_est(k-3), in_est(k-1), in_est(k-2), in_est(k-3)];
end

pLS3 = pinv(Phi3) * out_est;
G3 = tf(pLS3(4:6)', [1, pLS3(1:3)'], Tp);

%% IV - 3 rząd
x3 = lsim(G3, in_est, (0:N_est-1) * Tp);

Z3 = zeros(N_est, 6);
Z3(2, :) = [-x3(1), 0, 0, in_est(1), 0, 0];
Z3(3, :) = [-x3(2), -x3(1), 0, in_est(2), in_est(1), 0];

for k=4:N_est
    Z3(k, :) = [-x3(k-1), -x3(k-2), -x3(k-3), in_est(k-1), in_est(k-2), in_est(k-3)];
end

pIV3 = (Z3' * Phi3) \ Z3' * out_est;

G3IV = tf(pIV3(4:6)', [1, pIV3(1:3)'], Tp);

y_hat3 = zeros(N_wer, 1);
y_hat3(1) = out_wer(1);
y_hat3(2) = -pIV3(1) * out_wer(k-1) + pIV3(4) * in_wer(k-1);
y_hat3(3) = -pIV3(1) * out_wer(k-1) - pIV3(2) * out_wer(k-2) + pIV3(4) * in_wer(k-1) + pIV3(5) * in_wer(k-2);

for k=4:N_wer
    y_hat3(k) = -pIV3(1) * out_wer(k-1) - pIV3(2) * out_wer(k-2) - pIV3(3) * out_wer(k-3) + pIV3(4) * in_wer(k-1) + pIV3(5) * in_wer(k-2) + pIV3(6) * in_wer(k-3);
end

%% LS - 4 rząd
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) u(n-1) u(n-2) u(n-3) u(n-4)]
Phi4 = zeros(N_est, 8);
Phi4(2, :) = [-out_est(1), 0, 0, 0, in_est(1), 0, 0, 0];
Phi4(3, :) = [-out_est(2), -out_est(1), 0, 0, in_est(2), in_est(1), 0, 0];
Phi4(4, :) = [-out_est(3), -out_est(2), -out_est(1), 0, in_est(3), in_est(2), in_est(1), 0];

for k=5:N_est
    Phi4(k, :) = [-out_est(k-1), -out_est(k-2), -out_est(k-3), -out_est(k-4), in_est(k-1), in_est(k-2), in_est(k-3), in_est(k-4)];
end

pLS4 = pinv(Phi4) * out_est;
G4 = tf(pLS4(5:8)', [1, pLS4(1:4)'], Tp);

%% IV - 4 rząd
x4 = lsim(G4, in_est, (0:N_est-1) * Tp);

Z4 = zeros(N_est, 8);
Z4(2, :) = [-x4(1), 0, 0, 0, in_est(1), 0, 0, 0];
Z4(3, :) = [-x4(2), -x4(1), 0, 0, in_est(2), in_est(1), 0, 0];
Z4(4, :) = [-x4(3), -x4(2), -x4(1), 0, in_est(3), in_est(2), in_est(1), 0];

for k=5:N_est
    Z4(k, :) = [-x4(k-1), -x4(k-2), -x4(k-3), -x4(k-4), in_est(k-1), in_est(k-2), in_est(k-3), in_est(k-4)];
end

pIV4 = (Z4' * Phi4) \ Z4' * out_est;

G4IV = tf(pIV4(5:8)', [1, pIV4(1:4)'], Tp);

y_hat4 = zeros(N_wer, 1);
y_hat4(1) = out_wer(1);
y_hat4(2) = -pIV4(1) * out_wer(1) + pIV4(5) * in_wer(1);
y_hat4(3) = -pIV4(1) * out_wer(2) - pIV4(2) * out_wer(1) + pIV4(5) * in_wer(2) + pIV4(6) * in_wer(1);
y_hat4(4) = -pIV4(1) * out_wer(3) - pIV4(2) * out_wer(2) - pIV4(3) * out_wer(1) + pIV4(5) * in_wer(3) + pIV4(6) * in_wer(2) + pIV4(7) * in_wer(1);

for k=5:N_wer
    y_hat4(k) = -pIV4(1) * out_wer(k-1) - pIV4(2) * out_wer(k-2) - pIV4(3) * out_wer(k-3) - pIV4(4) * out_wer(k-4) + pIV4(5) * in_wer(k-1) + pIV4(6) * in_wer(k-2) + pIV4(7) * in_wer(k-3) + pIV4(8) * in_wer(k-4);
end

%% Porównanie odpowiedzi
t_wer = (0:N_wer-1)*Tp;
y3 = lsim(G3IV, in_wer, t_wer);
y4 = lsim(G4IV, in_wer, t_wer);

fig_arx_ls = figure('Position', [100 100 1200 900],...
                    'Name', 'Porównanie IV',...
                    'NumberTitle', 'off');
plot(t_wer, out_wer, '--');
hold on
plot(t_wer, y3);
plot(t_wer, y4);
plot(t_wer, y_hat3);
plot(t_wer, y_hat4);
legend('zmierzona odpowiedź',...
       'odp. modelu - 3 rząd',...
       'odp. modelu - 4 rząd',...
       'odp. predyktora - 3 rząd',...
       'odp. predyktora - 4 rząd');
xlabel('t [s]');
title('Porównanie IV');
grid on

%% J
J3 = mean((out_wer - y_hat3) .^ 2);
J4 = mean((out_wer - y_hat4) .^ 2);

%% Jfit
Jfit3 = (1 - norm(out_wer - y_hat3) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;
Jfit4 = (1 - norm(out_wer - y_hat4) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;

%% FPE
V3 = mean((out_est - Phi3 * pLS3) .^ 2);
FPE3 = V3 * (1 + 3 / N_est) / (1 - 3 / N_est);

V4 = mean((out_est - Phi4 * pLS4) .^ 2);
FPE4 = V4 * (1 + 4 / N_est) / (1 - 4 / N_est);

%% AIC
AIC3 = N_est * log(V3) + 2 * 3;
AIC4 = N_est * log(V4) + 2 * 4;

%% Współczynnik uwarunkowania macierzy MI

%TODO: Podmień wzory na IV
MI3 = 1/N_est * (Phi3' * Phi3);
cond3 = sqrt(max(eig(MI3' * MI3))) / sqrt(min(eig(MI3' * MI3)));

MI4 = 1/N_est * (Phi4' * Phi4);
cond4 = sqrt(max(eig(MI4' * MI4))) / sqrt(min(eig(MI4' * MI4)));

%% Podsumowanie
table([J3; J4],...
      [Jfit3; Jfit4],...
      [cond3; cond4],...
      [FPE3; FPE4],...
      [AIC3; AIC4],...
      'VariableNames', {'J', 'J_FIT', 'cond', 'FPE', 'AIC'},...
      'RowNames', {'M3', 'M4'})
