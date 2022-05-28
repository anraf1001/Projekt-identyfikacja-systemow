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

%% Identyfikacja - rząd 2
% phiT(n) = [-y(n-1) -y(n-2) u(n-1) u(n-2)]
Phi2 = zeros(N_est, 4);
Phi2(2, :) = [-out_est(1), 0, in_est(1), 0];

for k=3:N_est
    Phi2(k, :) = [-out_est(k-1), -out_est(k-2), in_est(k-1), in_est(k-2)];
end

pLS2 = pinv(Phi2) * out_est;
G2 = tf(pLS2(3:4)', [1, pLS2(1:2)'], Tp);

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
title('Porównanie LS');
grid on
