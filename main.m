close all; clear; clc;

%% Zainicjowanie zmiennych globalnych
Tp = 0.08;

%% Wczytanie oraz wyświetlenie danych pomiarowych
load('dane.mat');

fig_meas_data = figure('Position', [100 100 1200 900],...
                       'Name', 'Dane pomiarowe',...
                       'NumberTitle', 'off');

subplot(2, 1, 1);
plot((0:length(in) - 1) * Tp, in);
title('Sygnał wejściowy - napięcie nagrzewnicy');
xlabel('t [s]');
grid on

subplot(2, 1, 2);
plot((0:length(out) - 1) * Tp, out);
title('Sygnał wyjściowy - napięcie termopary');
xlabel('t [s]');
grid on
