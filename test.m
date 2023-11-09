

close all;
clear all;

%download variables to use
data0= readtable("C:\Users\malte\OneDrive\Desktop\thesis\Data\Clean\df_qoqld.csv", 'ReadVariableNames',true);

% switching series: output in data_SV
data_SV = data0(1:113,"ld_kpif_lag2");

Z=data_SV; % growth rate of inf
Z = table2array(Z);



% Create the updated time vector

time = 1:113


Th = prctile(Z, 60); % we do a 60/40 split
date = datetime(1995, 05, 31):calmonths(3):datetime(2023, 05, 31);

abo = Z > Th;

% Create a figure
figure;

% Plot the data points above and below the threshold
plot(date, Z, 'b');
hold on
xtickangle(45);
plot([date(1), date(end)], [Th Th], 'k--', 'LineWidth', 1); % threshold

% Find the indices where abo is equal to 1
above_indices = find(abo);

% Loop through the above_indices and create patches
% Loop through the above_indices and create patches
for i = 1:length(above_indices)
    if i < length(above_indices) && abo(above_indices(i))
        start_idx = above_indices(i);
        end_idx = above_indices(i+1);
        patch([date(start_idx), date(end_idx), date(end_idx), date(start_idx)], ...
              [0, 0, 10, 10], 'r');
    end
end

set(gca,'children',flipud(get(gca,'children')))
hold off;

%moving the threshold
Z0=Z-Th;


%% logistic regression

% Define the values for Z0

% Define the parameter theta (you can adjust this as needed)
theta = -20;  % Use a negative value for theta

% Calculate M_Z_t using the formula with a negative theta
M_Z_t = 1 ./ (1 + exp(theta * Z0));

[Z0_sorted, idx] = sort(Z0);
M_Z_t_sorted = M_Z_t(idx);

% Plot the sorted data with points and a line connecting them
plot(Z0_sorted, M_Z_t_sorted, 'b*');
hold on;
plot(Z0_sorted, M_Z_t_sorted, 'r-');
hold off;

xlabel('Z0');
ylabel('M_Z_t');
title('Plot of M_Z_t as a function of Z0 (Ordered)');
grid on;