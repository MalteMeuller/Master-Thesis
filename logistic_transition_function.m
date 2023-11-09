% Define the values for Z0
Z0 = linspace(-2.5, 2.5, 100);  % Adjust the range and number of points as needed

% Define the parameter theta (you can adjust this as needed)
theta = -12;  % Use a negative value for theta

% Calculate M_Z_t using the formula with a negative theta
M_Z_t = 1 ./ (1 + exp(theta * Z0));

% Plot the results
plot(Z0, M_Z_t, 'b', 'LineWidth', 2);
xlabel('Z0');
ylabel('M_Z_t');
title('Plot of M_Z_t as a function of Z0');
grid on;
%% 

% Define the values for Z0
Z0 = linspace(-2.5, 2.5, 100);  % Adjust the range and number of points as needed

% Define the parameter theta (you can adjust this as needed)
theta = -2;  % Use a negative value for theta

% Calculate M_Z_t using the formula with a negative theta
M_Z_t = 1 ./ (1 + exp(theta * Z0));

% Plot the results
plot(Z0, M_Z_t, 'b', 'LineWidth', 2);
xlabel('Z0');
ylabel('M_Z_t');
title('Plot of M_Z_t as a function of Z0');
grid on;


%%
%% logistic regression
clear all
close all
% Define the values for Z0
% Load the data
data0 = readtable("C:\thesis\Data\Clean\df_qoqld.csv", 'ReadVariableNames', true);

% Extract relevant data
Z0 = data0.ld_kpif_lag2;
stdz0 = std(Z0)
% Define the threshold
Th0 = prctile(Z0, 60); % 60/40 split

% Determine values above the threshold
abo0 = Z0 > Th0;

above = sum(abo0)
below = 113-above

Z1 = Z0-Th0
std = std(Z1)

% Define the parameter theta (you can adjust this as needed)
theta = -12;  % Use a negative value for theta

% Calculate M_Z_t using the formula with a negative theta
M_Z_t = 1 ./ (1 + exp(theta * Z1));
M_Z_t_nosc = 1 ./ (1 + exp(theta * Z1));

[Z1_sorted, idx] = sort(Z1);
M_Z_t_sorted = M_Z_t(idx);

% Plot the sorted data with points and a line connecting them
plot(Z1_sorted, M_Z_t_sorted, 'r-', 'LineWidth', 2);  % Increase line width
hold on;
plot(Z1_sorted, M_Z_t_sorted, 'k.', 'MarkerSize', 20);  % Increase marker size
hold off;
xlim([-1.6, 2.8]);

xlabel('CPIF_t_-_2');
ylabel('F(z_t)');
title('Logistic Smooth transition Function');
grid on;