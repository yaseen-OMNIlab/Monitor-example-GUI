clc; close all; clear;
%% load dataset
load corr_data.mat
%% plot
figure;plot(lags_stim,corr_data);
hold on
plot(lags_stim,mean(corr_data,2),'LineWidth',5);

%% plot shaded results
% Define the real time range
real_time_start = -5; % start of real time range in seconds
real_time_end = 5;    % end of real time range in seconds

% Calculate the scaling factor
scaling_factor = (5 - -5) / (max(lags_stim) - min(lags_stim));

% Calculate the offset
offset = real_time_start - min(lags_stim) * scaling_factor;

% Transform the time lags to real time
real_time = lags_stim * scaling_factor + offset;

% Calculate the mean and standard deviation across rows
mean_data = mean(corr_data, 2);
std_data = std(corr_data, 0, 2);

% Determine the sample size
n = size(data,1);

% Compute the standard error of the mean
sem_data = std_data / sqrt(n);

% Plot using shadedErrorBar
h = shadedErrorBar(real_time, mean_data, sem_data, 'lineProps', 'r');

% Find the minimum value in mean_data and its index
[min_mean, min_index] = min(mean_data);

% Get the corresponding standard deviation and time
min_sem = sem_data(min_index);
min_time = real_time(min_index);

% Display the results
fprintf('Minimum Mean Value: %f\n', min_mean);
fprintf('SEM at Minimum Mean: %f\n', min_sem);
fprintf('Time at Minimum Mean: %f seconds\n', min_time);

% Mark the minimum point on the plot
hold on;
plot(min_time, min_mean, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

% Customize the plot
xlabel('Time (seconds)');
ylabel('Mean Correlation');
title('Mean Correlation with SEM');
grid on;
legend('Mean ± SEM', 'Minimum Mean Value');