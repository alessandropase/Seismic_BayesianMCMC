%% Intro
clear;
close all;
clc;

%% Loading of simulations
load("test 1.mat");
mu_x_k1 = mu_x_k;
mu_y_k1 = mu_y_k;
load("test 2.mat");
mu_x_k2 = mu_x_k;
mu_y_k2 = mu_y_k;
load("test 3.mat");
mu_x_k3 = mu_x_k;
mu_y_k3 = mu_y_k;
load("test 3.mat");
mu_x_k4 = mu_x_k;
mu_y_k4 = mu_y_k;
load("test 5.mat");
mu_x_k5 = mu_x_k;
mu_y_k5 = mu_y_k;
load("test 6.mat");
mu_x_k6 = mu_x_k;
mu_y_k6 = mu_y_k;
load("test 6.mat");
mu_x_k7 = mu_x_k;
mu_y_k7 = mu_y_k;
load("test 8.mat");
mu_x_k8 = mu_x_k;
mu_y_k8 = mu_y_k;
load("test 8.mat");
mu_x_k9 = mu_x_k;
mu_y_k9 = mu_y_k;
load("test 10.mat");
mu_x_k10 = mu_x_k;
mu_y_k10 = mu_y_k;
load("test 11.mat");
mu_x_k11 = mu_x_k;
mu_y_k11 = mu_y_k;
load("test 12.mat");
mu_x_k12 = mu_x_k;
mu_y_k12 = mu_y_k;

k = 1:1:length(mu_x_k10);

figure;
hold on;
plot(k, mu_x_k1);
plot(k, mu_x_k2);
plot(k, mu_x_k3);
plot(k, mu_x_k4);
plot(k, mu_x_k5);
plot(k, mu_x_k6);
plot(k, mu_x_k7);
plot(k, mu_x_k8);
plot(k, mu_x_k9);
plot(k, mu_x_k10);
plot(k, mu_x_k11);
plot(k, mu_x_k12);

sigma_x = 72.03;
mu_x = 10.56;
r_x = @(n) sqrt(sigma_x^2./n) + mu_x;
r_x_min = @(n) -sqrt(sigma_x^2./n) + mu_x;
plot(k(2:end), r_x(k(2:end)), 'k', 'LineWidth', 2);
plot(k(2:end), r_x_min(k(2:end)), 'k', 'LineWidth', 2);
xlabel('Iteration', 'Interpreter','latex', 'FontSize', 15)
ylabel('$\mu_x$', 'Interpreter','latex', 'FontSize', 15)
grid on;
ylim([6 16]);
xlim([1, 50001])




figure;
hold on;
plot(k, mu_y_k1);
plot(k, mu_y_k2);
plot(k, mu_y_k3);
plot(k, mu_y_k4);
plot(k, mu_y_k5);
plot(k, mu_y_k6);
plot(k, mu_y_k7);
plot(k, mu_y_k8);
plot(k, mu_y_k9);
plot(k, mu_y_k10);
plot(k, mu_y_k11);
plot(k, mu_y_k12);
sigma_y = 108.64;
mu_y = 12.39;
r_y = @(n) sqrt(sigma_y^2./n) + mu_y;
r_y_min = @(n) -sqrt(sigma_y^2./n) + mu_y;
plot(k(2:end), r_y(k(2:end)), 'k', 'LineWidth', 2);
plot(k(2:end), r_y_min(k(2:end)), 'k', 'LineWidth', 2);
xlabel('Iteration', 'Interpreter','latex', 'FontSize', 15)
ylabel('$\mu_y$', 'Interpreter','latex', 'FontSize', 15)
grid on;
ylim([6 18]);
xlim([1, 50001]);



