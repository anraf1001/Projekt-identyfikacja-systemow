close all; clear; clc;

%% Zainicjowanie zmiennych globalnych
Tp = 0.08;

%% Wczytanie danych pomiarowych
load('dane.mat');
N = length(in);
t = (0:length(in) - 1) * Tp;

%% Podział danych pomiarowych na 2 podzbiory
in_est = in(1:fix(N/2));
in_wer = in(fix(N/2)+1:end);

out_est = out(1:fix(N/2));
out_wer = out(fix(N/2)+1:end);

N_est = length(in_est);
N_wer = length(in_wer);

%% Identyfikacja - rząd 1
% phiT(n) = [-y(n-1) u(n-1)]
Phi1 = zeros(N_est, 2);

for k=2:N_est
    Phi1(k, :) = [-out_est(k-1), in_est(k-1)];
end

pLS1 = pinv(Phi1) * out_est;
G1 = tf(pLS1(2), [1, pLS1(1)], Tp);

y_hat1 = zeros(N_wer, 1);
y_hat1(1) = out_wer(1);

for k=2:N_wer
    y_hat1(k) = -pLS1(1) * out_wer(k-1) + pLS1(2) * in_wer(k-1);
end

%% Identyfikacja - rząd 2
% phiT(n) = [-y(n-1) -y(n-2) u(n-1) u(n-2)]
Phi2 = zeros(N_est, 4);
Phi2(2, :) = [-out_est(1), 0, in_est(1), 0];

for k=3:N_est
    Phi2(k, :) = [-out_est(k-1), -out_est(k-2), in_est(k-1), in_est(k-2)];
end

pLS2 = pinv(Phi2) * out_est;
G2 = tf(pLS2(3:4)', [1, pLS2(1:2)'], Tp);

y_hat2 = zeros(N_wer, 1);
y_hat2(1) = out_wer(1);
y_hat2(2) = -pLS2(1) * out_wer(1) + pLS2(3) * in_wer(1);

for k=3:N_wer
    y_hat2(k) = -pLS2(1) * out_wer(k-1) - pLS2(2) * out_wer(k-2) + pLS2(3) * in_wer(k-1) + pLS2(4) * in_wer(k-2);
end

%% Identyfikacja - rząd 3
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) u(n-1) u(n-2) u(n-3)]
Phi3 = zeros(N_est, 6);
Phi3(2, :) = [-out_est(1), 0, 0, in_est(1), 0, 0];
Phi3(3, :) = [-out_est(2), -out_est(1), 0, in_est(2), in_est(1), 0];

for k=4:N_est
    Phi3(k, :) = [-out_est(k-1), -out_est(k-2), -out_est(k-3), in_est(k-1), in_est(k-2), in_est(k-3)];
end

pLS3 = pinv(Phi3) * out_est;
G3 = tf(pLS3(4:6)', [1, pLS3(1:3)'], Tp);

y_hat3 = zeros(N_wer, 1);
y_hat3(1) = out_wer(1);
y_hat3(2) = -pLS3(1) * out_wer(1) + pLS3(4) * in_wer(1);
y_hat3(3) = -pLS3(1) * out_wer(2) - pLS3(2) * out_wer(1) + pLS3(4) * in_wer(2) + pLS3(5) * in_wer(1);

for k=4:N_wer
    y_hat3(k) = -pLS3(1) * out_wer(k-1) - pLS3(2) * out_wer(k-2) - pLS3(3) * out_wer(k-3) + pLS3(4) * in_wer(k-1) + pLS3(5) * in_wer(k-2) + pLS3(6) * in_wer(k-3);
end

%% Identyfikacja - rząd 4
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

y_hat4 = zeros(N_wer, 1);
y_hat4(1) = out_wer(1);
y_hat4(2) = -pLS4(1) * out_wer(1) + pLS4(5) * in_wer(1);
y_hat4(3) = -pLS4(1) * out_wer(2) - pLS4(2) * out_wer(1) + pLS4(5) * in_wer(2) + pLS4(6) * in_wer(1);
y_hat4(4) = -pLS4(1) * out_wer(3) - pLS4(2) * out_wer(2) - pLS4(3) * out_wer(1) + pLS4(5) * in_wer(3) + pLS4(6) * in_wer(2) + pLS4(7) * in_wer(1);

for k=5:N_wer
    y_hat4(k) = -pLS4(1) * out_wer(k-1) - pLS4(2) * out_wer(k-2) - pLS4(3) * out_wer(k-3) - pLS4(4) * out_wer(k-4) + pLS4(5) * in_wer(k-1) + pLS4(6) * in_wer(k-2) + pLS4(7) * in_wer(k-3) + pLS4(8) * in_wer(k-4);
end

%% Identyfikacja - rząd 5
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) -y(n-5) u(n-1) u(n-2) u(n-3) u(n-4) u(n-5)]
Phi5 = zeros(N_est, 10);
Phi5(2, :) = [-out_est(1), 0, 0, 0, 0, in_est(1), 0, 0, 0, 0];
Phi5(3, :) = [-out_est(2), -out_est(1), 0, 0, 0, in_est(2), in_est(1), 0, 0, 0];
Phi5(4, :) = [-out_est(3), -out_est(2), -out_est(1), 0, 0, in_est(3), in_est(2), in_est(1), 0, 0];
Phi5(5, :) = [-out_est(4), -out_est(3), -out_est(2), -out_est(1), 0, in_est(4), in_est(3), in_est(2), in_est(1), 0];

for k=6:N_est
    Phi5(k, :) = [-out_est(k-1), -out_est(k-2), -out_est(k-3), -out_est(k-4), -out_est(k-5), in_est(k-1), in_est(k-2), in_est(k-3), in_est(k-4), in_est(k-5)];
end

pLS5 = pinv(Phi5) * out_est;
G5 = tf(pLS5(6:10)', [1, pLS5(1:5)'], Tp);

y_hat5 = zeros(N_wer, 1);
y_hat5(1) = out_wer(1);
y_hat5(2) = -pLS5(1) * out_wer(1) + pLS5(6) * in_wer(1);
y_hat5(3) = -pLS5(1) * out_wer(2) - pLS5(2) * out_wer(1) + pLS5(6) * in_wer(2) + pLS5(7) * in_wer(1);
y_hat5(4) = -pLS5(1) * out_wer(3) - pLS5(2) * out_wer(2) - pLS5(3) * out_wer(1) + pLS5(6) * in_wer(3) + pLS5(7) * in_wer(2) + pLS5(8) * in_wer(1);
y_hat5(5) = -pLS5(1) * out_wer(4) - pLS5(2) * out_wer(3) - pLS5(3) * out_wer(2) - pLS5(4) * out_wer(1) + pLS5(6) * in_wer(4) + pLS5(7) * in_wer(3) + pLS5(8) * in_wer(2) + pLS5(9) * in_wer(1);

for k=6:N_wer
    y_hat5(k) = -pLS5(1) * out_wer(k-1) - pLS5(2) * out_wer(k-2) - pLS5(3) * out_wer(k-3) - pLS5(4) * out_wer(k-4) - pLS5(5) * out_wer(k-5) + pLS5(6) * in_wer(k-1) + pLS5(7) * in_wer(k-2) + pLS5(8) * in_wer(k-3) + pLS5(9) * in_wer(k-4) + pLS5(10) * in_wer(k-5);
end

%% Porównanie odpowiedzi
t_wer = (0:N_wer-1)*Tp;
y1 = lsim(G1, in_wer, t_wer);
y2 = lsim(G2, in_wer, t_wer);
y3 = lsim(G3, in_wer, t_wer);
y4 = lsim(G4, in_wer, t_wer);
y5 = lsim(G5, in_wer, t_wer);

fig_arx_ls = figure('Position', [100 100 1200 900],...
                    'Name', 'Porównanie LS',...
                    'NumberTitle', 'off');
plot(t_wer, out_wer, '--');
hold on
plot(t_wer, y1);
plot(t_wer, y2);
plot(t_wer, y3);
plot(t_wer, y4);
plot(t_wer, y5);
legend('zmierzona odpowiedź',...
       'odp. modelu - 1 rząd',...
       'odp. modelu - 2 rząd',...
       'odp. modelu - 3 rząd',...
       'odp. modelu - 4 rząd',...
       'odp. modelu - 5 rząd');
xlabel('t [s]');
title('Porównanie LS - odpowiedź modelu');
grid on

fig_arx_ls_pred = figure('Position', [100 100 1200 900],...
                         'Name', 'Porównanie LS',...
                         'NumberTitle', 'off');
plot(t_wer, out_wer, '--');
hold on
plot(t_wer, y_hat1);
plot(t_wer, y_hat2);
plot(t_wer, y_hat3);
plot(t_wer, y_hat4);
plot(t_wer, y_hat5);
legend('zmierzona odpowiedź',...
       'odp. predyktora - 1 rząd',...
       'odp. predyktora - 2 rząd',...
       'odp. predyktora - 3 rząd',...
       'odp. predyktora - 4 rząd',...
       'odp. predyktora - 5 rząd');
xlabel('t [s]');
title('Porównanie LS - predyktor jednokrokowy');
grid on

%% J
J1 = mean((out_wer - y_hat1) .^ 2);
J2 = mean((out_wer - y_hat2) .^ 2);
J3 = mean((out_wer - y_hat3) .^ 2);
J4 = mean((out_wer - y_hat4) .^ 2);
J5 = mean((out_wer - y_hat5) .^ 2);

%% Jfit
Jfit1 = (1 - norm(out_wer - y_hat1) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;
Jfit2 = (1 - norm(out_wer - y_hat2) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;
Jfit3 = (1 - norm(out_wer - y_hat3) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;
Jfit4 = (1 - norm(out_wer - y_hat4) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;
Jfit5 = (1 - norm(out_wer - y_hat5) / norm(out_wer - mean(out_wer)*ones(size(out_wer)))) * 100;

%% FPE
V1 = mean((out_est - Phi1 * pLS1) .^ 2);
FPE1 = V1 * (1 + 1 / N_est) / (1 - 1 / N_est);

V2 = mean((out_est - Phi2 * pLS2) .^ 2);
FPE2 = V2 * (1 + 2 / N_est) / (1 - 2 / N_est);

V3 = mean((out_est - Phi3 * pLS3) .^ 2);
FPE3 = V3 * (1 + 3 / N_est) / (1 - 3 / N_est);

V4 = mean((out_est - Phi4 * pLS4) .^ 2);
FPE4 = V4 * (1 + 4 / N_est) / (1 - 4 / N_est);

V5 = mean((out_est - Phi5 * pLS5) .^ 2);
FPE5 = V5 * (1 + 5 / N_est) / (1 - 5 / N_est);

%% AIC

AIC1 = N_est * log(V1) + 2 * 1;
AIC2 = N_est * log(V2) + 2 * 2;
AIC3 = N_est * log(V3) + 2 * 3;
AIC4 = N_est * log(V4) + 2 * 4;
AIC5 = N_est * log(V5) + 2 * 5;

%% Współczynnik uwarunkowania macierzy MI
MI1 = 1/N_est * (Phi1' * Phi1);
cond1 = sqrt(max(eig(MI1' * MI1))) / sqrt(min(eig(MI1' * MI1)));

MI2 = 1/N_est * (Phi2' * Phi2);
cond2 = sqrt(max(eig(MI2' * MI2))) / sqrt(min(eig(MI2' * MI2)));

MI3 = 1/N_est * (Phi3' * Phi3);
cond3 = sqrt(max(eig(MI3' * MI3))) / sqrt(min(eig(MI3' * MI3)));

MI4 = 1/N_est * (Phi4' * Phi4);
cond4 = sqrt(max(eig(MI4' * MI4))) / sqrt(min(eig(MI4' * MI4)));

MI5 = 1/N_est * (Phi5' * Phi5);
cond5 = sqrt(max(eig(MI5' * MI5))) / sqrt(min(eig(MI5' * MI5)));

%% Podsumowanie
table([J1; J2; J3; J4; J5],...
      [Jfit1; Jfit2; Jfit3; Jfit4; Jfit5],...
      [cond1; cond2; cond3; cond4; cond5],...
      [FPE1; FPE2; FPE3; FPE4; FPE5],...
      [AIC1; AIC2; AIC3; AIC4; AIC5],...
      'VariableNames', {'J', 'J_FIT', 'cond', 'FPE', 'AIC'},...
      'RowNames', {'M1', 'M2', 'M3', 'M4', 'M5'})
