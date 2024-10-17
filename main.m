%% intro
clc;
clear;
close all;

%% Datas
sens_pos = [3, 15;
            3, 16;
            4, 15;
            4, 16;
            5, 15;
            5, 16];                                                         %  [ km ]
t_obs = [3.12, 3.26, 2.98, 3.12, 2.84, 2.98];                               %  [ s ]

sigma = 0.1;                                                                %  [ s ]
v = 5;                                                                      % [ km/s ]

%% Point 1.
x_j = sens_pos(:, 1);                                                       % Vector of x position of the sensors
y_j = sens_pos(:, 2);                                                       % Vector of y position of the sensors
n_sens = length(x_j);                                                       % Number of the sensors

% Plot variables
n = 500;
xx = linspace(0, 20, n);
yy = linspace(0, 20, n);
[X, Y] = meshgrid(xx, yy);

% Posterior probability density function computation
pi_prior = 1;                                                               % Prior PDF
pi_post = @(x, y) pi_prior * likelihood(x, y, x_j, y_j, t_obs, v, sigma);   % Posterior PDF - Bayesian formulation

% Normalization of the posterior
norm_const = integral2(@(x, y) pi_post(x, y), -50, 50, -50, 50);
pi_post = @(x,y) pi_post(x,y)/norm_const;                                   % The constant of the uniform prior PDF is chosen s.t. the integral in the domain of the posterior is 1

% Evaluating the posterior PDF
PP = zeros(length(xx));
for a = 1:length(xx)
    for j = 1:length(yy)
        PP(j, a) = pi_post(xx(a), yy(j));
    end
end

%% Creating the plot of the posterior
figure
s = surf(X, Y, PP);                                                         % Surface plot
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)
s.EdgeColor = 'none';
colormap jet

figure
contour(X, Y, PP, 300, 'ShowText','off')                                    % Contour plot
hold on;
plot(x_j, y_j, 'o', 'LineWidth', 3)
grid on;
legend('Contour lines of posterior PDF', 'Sensor positions', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10)
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)

%% Statistical moments 
mu_x = integral2(@(x, y) x .* pi_post(x, y), -50, 50, -50, 50);             % mu_x
mu_y = integral2(@(x, y) y .* pi_post(x, y), -50, 50, -50, 50);             % mu_y
var_x = integral2(@(x, y) (x-mu_x).^2 .* pi_post(x, y), -50, 50, -50, 50);  % var_x
var_y = integral2(@(x, y) (y-mu_y).^2 .* pi_post(x, y), -50, 50, -50, 50);  % var_y

fprintf('Mean for x: %.2f\n', mu_x);
fprintf('Mean for y: %.2f\n', mu_y);
fprintf('Variance for x: %.2f\n', var_x);
fprintf('Variance for y: %.2f\n', var_y);

%% 3. Sampling
% Choice of starting point: the starting point is chosen randomly in the domain
% [0,20] x [0, 20] such that is posterior is different from zero
theta_0 = [20*rand(1), 20*rand(1)];
toll = 1e-4;
while(pi_post(theta_0(1), theta_0(2))<toll)
    theta_0 = [20*rand(1), 20*rand(1)];
end

k_max = 100000;                                                             % Number of considered samples
theta = zeros(2, k_max+1);                                                  % Initializing the samples vector
theta(:, 1) = theta_0;
mu_x_k = zeros(1, k_max/2+1);                                               % Initializing the mean of x vector
mu_y_k = zeros(1, k_max/2+1);                                               % Initializing the mean of y vector

scaling = 1;
a = 1;                                                                      % Counter for the accepted samples
b = 0;
for k = 2 : k_max+1
    theta_old = theta(:, k-1);
    SIGMA = scaling*[var_x, 0;
             0, var_y];                                                     % The matrix is chosen as a diagonal matrix with the values of the variance
    m = 2;                                                                  % Dimension of the problem
    C = 2.4^2 / m * SIGMA;                                                  % Gelman et al. proposed [C] matrix
    theta_hat = mvnrnd(theta_old, C, 1);                                    % Sampling from a multivariate gaussian distribution
    alpha = min(1, pi_post(theta_hat(1), theta_hat(2)) / pi_post(theta_old(1), theta_old(2))); % Acceptance probability
    u = rand(1);
    if alpha >= u                                                           % New sample is accepted
        theta(:, k) = theta_hat;
        a = a + 1;
    else                                                                    % New sample is discarded
        theta(:, k) = theta_old;
    end
    
    if k > k_max/2                                                          % The mean is calculated according to the burn-in strategy
        b = b+1;
        mu_x_k(b) = 1/b*sum(theta(1, k_max/2+1 : k));
        mu_y_k(b) = 1/b*sum(theta(2, k_max/2+1 : k));
    end
end

acc_rateo = a/k;                                                            % Acceptance rateo
fprintf('Acceptance rateo: %.2f\n', acc_rateo);

% Plot of the mean value wrt iterations (Burn-in)
figure;
plot(k_max/2+1:k_max/2+b, mu_x_k, 'LineWidth',2);
hold on;
grid on;
plot(k_max/2+1:k_max/2+b, mu_x*ones(size(mu_x_k)), '--', 'LineWidth',2);
xlim([k_max/2+1, k_max/2+b+2])
xlabel('Iterations', 'FontSize', 15, 'Interpreter','latex')
ylabel('Mean x value', 'FontSize',15, 'Interpreter','latex')
legend('Value obtained from MCMC', 'Numerically computed value from the Bayesian posterior', 'Location','best')

figure;
plot(k_max/2+1:k_max/2+b, mu_y_k, 'LineWidth',2);
hold on;
grid on;
plot(k_max/2+1:k_max/2+b, mu_y*ones(size(mu_y_k)), '--', 'LineWidth',2);
xlim([k_max/2+1, k_max/2+b+2])
xlabel('Iterations', 'FontSize', 15, 'Interpreter','latex')
ylabel('Mean y value', 'FontSize',15, 'Interpreter','latex')
legend('Value obtained from MCMC', 'Numerically computed value from the Bayesian posterior', 'Location','best', 'Interpreter', 'latex', 'FontSize', 10)


% Plot of the Markov-chain on top of the contour plot of the posterior
figure;
contour(X, Y, PP, 30, 'LineWidth', 1.5);
hold on;
plot(theta(1, :), theta(2, :), 'o-')
grid on;
legend('Contour lines of posterior PDF', 'Markov-chain trajectory', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10)
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)

% Plot of a 3D histogram for the samples
figure;
hist3(theta','Nbins',[200 200], 'CDataMode','auto','FaceColor','interp', 'Edgecolor', 'none');
colormap jet
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)

%% 4. Convergence study of the MCMC
N = 10;                                                                     % Number of MCMC considered
mean_x = zeros(1,N);                                                        % Vector pre-allocation
mean_y = zeros(1,N);                                                        % Vector pre-allocation
variance_x = zeros(1,N);
variance_y = zeros(1,N);
for j=1:N                                                                   % For N times the sampling is performed as before
    theta_0 = [20*rand(1), 20*rand(1)];
    toll = 1e-4;
    while(pi_post(theta_0(1), theta_0(2))<toll)
        theta_0 = [20*rand(1), 20*rand(1)];
    end
    k_max = 1000;                                                           % Number of samples in a MCMC
    theta = zeros(2, k_max+1);
    theta(:, 1) = theta_0;
    i = 1;
    for k = 2 : k_max+1
        theta_old = theta(:, k-1);
        SIGMA = 0.001*[var_x 0;
                       0 var_y];                                            % Variance calculated numerically
        m = 2;
        C = 2.4^2 / m * SIGMA;
        theta_hat = mvnrnd(theta_old, C, 1);
        alpha = min(1, pi_post(theta_hat(1), theta_hat(2)) / pi_post(theta_old(1), theta_old(2)));
        u = rand(1);
        if alpha >= u
            theta(:, k) = theta_hat;
            i = i + 1;
        else
            theta(:, k) = theta_old;
        end
    end
    k_new = k_max/2 + 1;
    theta = theta(:,k_new:end);
    mean_x(j) = mean(theta(1,:));                                           
    mean_y(j) = mean(theta(2,:));
    variance_x(j) = var(theta(1,:));
    variance_y(j) = var(theta(2,:));
end                                                                         % The sampling is the same as Section 3. but the means and the variances are computed
B_x = var(mean_x);                                                          % Variance of the means along x
B_y = var(mean_y);                                                          % Variance of the means along y
n = length(theta);                                                          
sigma_x_hat_sqrd = (n-1).*variance_x./n+B_x/n;                              
sigma_y_hat_sqrd = (n-1).*variance_y./n+B_y/n;
mu_x_hat = mean(mean_x);
mu_y_hat = mean(mean_y);
V_x_hat = sigma_x_hat_sqrd + B_x/N/n;
V_y_hat = sigma_y_hat_sqrd + B_y/N/n;
d_x = 2.*V_x_hat.^2/var(V_x_hat);
d_y = 2.*V_y_hat.^2/var(V_y_hat);

% G-R analysis
R_hat = zeros(2,N);                                                         % Vector pre-allocation
for i=1:N
    R_hat(1, i) = sqrt(((d_x(i)+3)*V_x_hat(i))./((d_x(i)+1)*variance_x(i)));
    R_hat(2, i) = sqrt(((d_y(i)+3)*V_y_hat(i))./((d_y(i)+1)*variance_y(i)));
end
R_hat_comp = max(max(R_hat));

%% 5. Sampling with adaptive scaling
% Choice of starting point: the starting point is chosen randomly in the domain
% [0,20] x [0, 20] such that is posterior is different from zero
theta_0 = [20*rand(1), 20*rand(1)];
toll = 1e-4;
while(pi_post(theta_0(1), theta_0(2))<toll)
    theta_0 = [20*rand(1), 20*rand(1)];
end
close all;
k_max = 200000;                                                             % Number of considered samples
theta = zeros(2, k_max+1);                                                  % Initializing the samples vector
theta(:, 1) = theta_0;
mu_x_k = zeros(1, k_max/2+1);                                               % Initializing the mean of x vector
mu_y_k = zeros(1, k_max/2+1);                                               % Initializing the mean of y vector
k_trust = 100;
scaling = 0.1;                                                              % Initial value fo the scaling
acc_opt = 0.23;
a = 1;                                                                      % Counter for the accepted samples
b = 0;
for k = 2 : k_max+1
    theta_old = theta(:, k-1);
    SIGMA = scaling*[var_x, 0;
             0, var_y];                                                     % The matrix is chosen as a diagonal matrix with the values of the variance
    m = 2;                                                                  % Dimension of the problem
    C = 2.4^2 / m * SIGMA;                                                  % Gelman et al. proposed [C] matrix
    theta_hat = mvnrnd(theta_old, C, 1);                                    % Sampling from a multivariate gaussian distribution
    alpha = min(1, pi_post(theta_hat(1), theta_hat(2)) / pi_post(theta_old(1), theta_old(2))); % Acceptance probability
    u = rand(1);
    if alpha >= u                                                           % New sample is accepted
        theta(:, k) = theta_hat;
        a = a + 1;
    else                                                                    % New sample is discarded
        theta(:, k) = theta_old;
    end
    
    if k > k_max/2                                                          % The mean is calculated according to the burn-in strategy
        b = b+1;
        mu_x_k(b) = 1/b*sum(theta(1, k_max/2+1 : k));
        mu_y_k(b) = 1/b*sum(theta(2, k_max/2+1 : k));
    end

    acceptance = a / k;

    if k >= k_trust && acceptance < acc_opt
        scaling = scaling *(1-abs(acc_opt-acceptance));
    elseif k >= k_trust && acceptance > acc_opt
        scaling = scaling *(1+abs(acc_opt-acceptance));
    end
    if scaling > 0.1
        scaling = 0.1;
    end

end

acc_rateo = a/k;                                                            % Acceptance rateo
fprintf('Acceptance rateo: %.2f\n', acc_rateo);

% Plot of the mean value wrt iterations (Burn-in)
figure;
plot(k_max/2+1:k_max/2+b, mu_x_k, 'LineWidth',2);
hold on;
grid on;
plot(k_max/2+1:k_max/2+b, mu_x*ones(size(mu_x_k)), '--', 'LineWidth',2);
xlim([k_max/2+1, k_max/2+b+2])
xlabel('Iterations', 'FontSize', 15, 'Interpreter','latex')
ylabel('Mean x value', 'FontSize',15, 'Interpreter','latex')
legend('Value obtained from MCMC', 'Numerically computed value from the Bayesian posterior', 'Location','best')

figure;
plot(k_max/2+1:k_max/2+b, mu_y_k, 'LineWidth',2);
hold on;
grid on;
plot(k_max/2+1:k_max/2+b, mu_y*ones(size(mu_y_k)), '--', 'LineWidth',2);
xlim([k_max/2+1, k_max/2+b+2])
xlabel('Iterations', 'FontSize', 15, 'Interpreter','latex')
ylabel('Mean y value', 'FontSize',15, 'Interpreter','latex')
legend('Value obtained from MCMC', 'Numerically computed value from the Bayesian posterior', 'Location','best', 'Interpreter', 'latex', 'FontSize', 10)

% Plot of the Markov-chain on top of the contour plot of the posterior
figure;
contour(X, Y, PP, 30, 'LineWidth', 1.5);
hold on;
plot(theta(1, :), theta(2, :), 'o-')
grid on;
legend('Contour lines of posterior PDF', 'Markov-chain trajectory', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10)
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)

% Plot of a 3D histogram for the samples
figure;
hist3(theta','Nbins',[200 200], 'CDataMode','auto','FaceColor','interp', 'Edgecolor', 'none');
colormap jet
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)

%% 6. Sampling without knowing sigma + adaptive coefficient
% Choice of starting point: the starting point is chosen randomly in the domain
% [0,20] x [0, 20] such that is posterior is different from zero
theta_0 = [20*rand(1), 20*rand(1)];
toll = 1e-4;
while(pi_post(theta_0(1), theta_0(2))<toll)
    theta_0 = [20*rand(1), 20*rand(1)];
end

k_max = 100000;                                                             % Number of considered samples
theta = zeros(2, k_max+1);                                                  % Initializing the samples vector
theta(:, 1) = theta_0;
mu_x_k = zeros(1, k_max/2+1);                                               % Initializing the mean of x vector
mu_y_k = zeros(1, k_max/2+1);                                               % Initializing the mean of y vector

a = 1;                                                                      % Counter for the accepted samples
b = 0;

C_0 = eye(2);                                                               % Initial guess of the matrix ==> I
%C_0 = [var_y, 0;
%           0, var_x];                                                      % UNCOMMENT to use a very wrong initial guess of the matrix 
k_trust = 100;                                                              % Length of the non adaptive interval
eps = 0;
acc_opt = 0.23;                                                             % Optimal acceptance rateo
scaling = 0.1;                                                              % Starting scaling value

for k = 2 : k_max+1
    theta_old = theta(:, k-1);
    if k < k_trust
        SIGMA = C_0;
        m = 2;
        C = 2.4^2 / m * SIGMA;
    else
        m = 2;
        C = scaling * 2.4^2 / m * cov(theta(:, 1:k-1)') + eps*eye(2);       % Adaptive value of the matrix only after the fixed iterations
    end
   
    theta_hat = mvnrnd(theta_old, C, 1);
    alpha = min(1, pi_post(theta_hat(1), theta_hat(2)) / pi_post(theta_old(1), theta_old(2)));
    u = rand(1);
    if alpha >= u
        theta(:, k) = theta_hat;
        a = a + 1;
    else
        theta(:, k) = theta_old;
    end
    if k > k_max/2                                                          % The mean is calculated according to the burn-in strategy
        b = b+1;
        mu_x_k(b) = 1/b*sum(theta(1, k_max/2+1 : k));
        mu_y_k(b) = 1/b*sum(theta(2, k_max/2+1 : k));
    end
    acceptance = a / k;

    if k >= k_trust && acceptance < acc_opt
        scaling = scaling *(1-abs(acc_opt-acceptance));
    elseif k >= k_trust && acceptance > acc_opt
        scaling = scaling *(1+abs(acc_opt-acceptance));
    end
    if scaling > 0.5
        scaling = 0.5;
    end
end

acc_rateo = a/k;                                                            % Acceptance rateo
fprintf('Acceptance rateo: %.2f\n', acc_rateo);

% Plot of the mean value wrt iterations (Burn-in)
figure;
plot(k_max/2+1:k_max/2+b, mu_x_k, 'LineWidth',2);
hold on;
grid on;
plot(k_max/2+1:k_max/2+b, mu_x*ones(size(mu_x_k)), '--', 'LineWidth',2);
xlim([k_max/2+1, k_max/2+b+2])
xlabel('Iterations', 'FontSize', 15, 'Interpreter','latex')
ylabel('Mean x value', 'FontSize',15, 'Interpreter','latex')
legend('Value obtained from MCMC', 'Numerically computed value from the Bayesian posterior', 'Location','best')

figure;
plot(k_max/2+1:k_max/2+b, mu_y_k, 'LineWidth',2);
hold on;
grid on;
plot(k_max/2+1:k_max/2+b, mu_y*ones(size(mu_y_k)), '--', 'LineWidth',2);
xlim([k_max/2+1, k_max/2+b+2])
xlabel('Iterations', 'FontSize', 15, 'Interpreter','latex')
ylabel('Mean y value', 'FontSize',15, 'Interpreter','latex')
legend('Value obtained from MCMC', 'Numerically computed value from the Bayesian posterior', 'Location','best', 'Interpreter', 'latex', 'FontSize', 10)

% Plot of the Markov-chain on top of the contour plot of the posterior
figure;
contour(X, Y, PP, 30, 'LineWidth', 1.5);
hold on;
plot(theta(1, :), theta(2, :), 'o-')
grid on;
legend('Contour lines of posterior PDF', 'Markov-chain trajectory', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10)
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)

% Plot of a 3D histogram for the samples
figure;
hist3(theta','Nbins',[200 200], 'CDataMode','auto','FaceColor','interp', 'Edgecolor', 'none');
colormap jet
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('y', 'Interpreter','latex', 'FontSize',15)