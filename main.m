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
