% =========================================================================
% @file    neuron_weight_prediction_roc_auc.m
% @brief   1000-fold train-test validation & ROC/AUC analysis for 40-area neuron weight prediction
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key parameters: Distance threshold=33.5, Gradient threshold=0.19 (40x40 version)
% =========================================================================

%% ===================== Part 1: Data Loading & Preprocessing =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;
rng('shuffle');  % Set random seed for reproducible train-test splits

% Load core datasets (critical dependencies)

load('MatF_SLN2021.mat');               % 91x91 neuron weight matrix (ground truth)
load('lp19_gd_fitnew.mat');             % Gradient fit data (gdfit10: 91x1 gradient values)
load('Distance_Matrix.mat');            % Pairwise distance matrix (dist_m: 91x91)
load('Normal40area_40_40.mat');         % Normalization parameters (nomal: 4x4 matrix)

% Preprocess weight matrix: retain only non-zero columns (40 valid areas)
weight_matrix_91x91 = MatF;
non_zero_columns = any(weight_matrix_91x91);  % Logical mask for non-zero columns
valid_area_indices = find(non_zero_columns);  % Indices of 40 valid brain areas

% Create 91x91 matrix with only 40x40 valid area submatrix (others zero)
weight_matrix_40area = zeros(91, 91);
weight_matrix_40area(valid_area_indices, valid_area_indices) = weight_matrix_91x91(valid_area_indices, valid_area_indices);

% Build full edge matrix (all 91x91 non-self area pairs: 8190 total edges)
clear full_edge_matrix;
edge_count = 0;

for source_area = 1:91
    for target_area = 1:91
        if source_area ~= target_area  % Exclude self-connections (i==k)
            edge_count = edge_count + 1;
            full_edge_matrix(edge_count, 1) = source_area;                          % Source area index
            full_edge_matrix(edge_count, 2) = target_area;                          % Target area index
            full_edge_matrix(edge_count, 3) = abs(log(gdfit10(source_area,1)/gdfit10(target_area,1)));  % Gradient ratio
            full_edge_matrix(edge_count, 4) = dist_m(source_area, target_area);      % Pairwise distance
        end
    end
end


% Define 4 bins (split by distance + gradient thresholds for 40x40 analysis)
distance_threshold = 33.5;  % Distance cutoff (mm)
gradient_threshold = 0.19;  % Gradient ratio cutoff

% Assign full edges to 4 bins (consistent thresholds for all iterations)
bin_indices_full{1} = find((full_edge_matrix(:,4) <= distance_threshold) & (full_edge_matrix(:,3) <= gradient_threshold));
bin_indices_full{2} = find((full_edge_matrix(:,4) > distance_threshold) & (full_edge_matrix(:,3) <= gradient_threshold));
bin_indices_full{3} = find((full_edge_matrix(:,4) <= distance_threshold) & (full_edge_matrix(:,3) > gradient_threshold));
bin_indices_full{4} = find((full_edge_matrix(:,4) > distance_threshold) & (full_edge_matrix(:,3) > gradient_threshold));

% Preprocess full edge bins with global normalization (for prediction)
clear full_bin_data;
for bin_idx = 1:4
    full_bin_data{bin_idx} = full_edge_matrix(bin_indices_full{bin_idx}, :);
    
    % Apply global normalization (from Normal40area_40_40.mat)
    full_bin_data{bin_idx}(:,8) = (full_bin_data{bin_idx}(:,3) - nomal(bin_idx,1)) / nomal(bin_idx,2);  % Normalized gradient
    full_bin_data{bin_idx}(:,9) = (full_bin_data{bin_idx}(:,4) - nomal(bin_idx,3)) / nomal(bin_idx,4);  % Normalized distance
end
fprintf('4 bins normalized: Bin 1=%d, Bin 2=%d, Bin 3=%d, Bin 4=%d edges\n', ...
    size(bin_indices_full{1},1), size(bin_indices_full{2},1), size(bin_indices_full{3},1), size(bin_indices_full{4},1));

%% ===================== Part 2: 1000 Iterations of Train-Test Split =====================
% Initialize storage variables
num_valid_iterations = 0;    % Counter for valid train-test splits
max_attempts = 5000;         % Max attempts to find valid splits (avoid infinite loop)
total_attempts = 0;          % Total attempts (valid + invalid)
min_bin_samples = 5;         % Minimum samples per bin for valid regression

% Cell arrays to store results (1000 iterations × metrics)
train_test_matrices = cell(1000, 4);  % [train_matrix, test_matrix, train_areas, test_areas]
regression_models = cell(1000, 4);    % Linear regression models (4 bins per iteration)
predicted_log_weights = cell(1000, 5);% Predicted log weights (4 bins + combined)

% Main loop: 1000 valid train-test splits (50/50 split of 40 areas)

while num_valid_iterations < 1000 && total_attempts < max_attempts
    total_attempts = total_attempts + 1;
    
    % Randomly split 40 valid areas into 2 train/test subsets (20 each)
    random_area_idx = randperm(40);
    train_area_indices = valid_area_indices(random_area_idx(1:20));
    test_area_indices = valid_area_indices(random_area_idx(21:40));
    
    % Build training matrix (zero out train areas to isolate test data)
    train_matrix = weight_matrix_40area;
    train_matrix(:, train_area_indices) = 0;  % Zero out train columns
    train_matrix(train_area_indices, :) = 0;  % Zero out train rows
    
    % Extract training data (non-zero entries from train matrix)
    clear train_feature_matrix;
    train_row_count = 0;
    
    for source_area = 1:91
        for target_area = 1:91
            if train_matrix(source_area, target_area) ~= 0
                train_row_count = train_row_count + 1;
                train_feature_matrix(train_row_count, 1) = source_area;
                train_feature_matrix(train_row_count, 2) = target_area;
                train_feature_matrix(train_row_count, 3) = log(train_matrix(source_area, target_area));  % Log-transformed weight
                train_feature_matrix(train_row_count, 4) = train_matrix(source_area, target_area);      % Raw weight
                train_feature_matrix(train_row_count, 5) = abs(log(gdfit10(source_area,1)/gdfit10(target_area,1)));  % Gradient
                train_feature_matrix(train_row_count, 6) = dist_m(source_area, target_area);            % Distance
            end
        end
    end
    
    % Split training data into 4 bins (same thresholds as full data)
    ind_dist_low = find(train_feature_matrix(:,6) <= distance_threshold);
    ind_dist_high = find(train_feature_matrix(:,6) > distance_threshold);
    train_dist_low = train_feature_matrix(ind_dist_low, :);
    train_dist_high = train_feature_matrix(ind_dist_high, :);
    
    ind_grad_low_low = find(train_dist_low(:,5) <= gradient_threshold);
    ind_grad_high_low = find(train_dist_low(:,5) > gradient_threshold);
    ind_grad_low_high = find(train_dist_high(:,5) <= gradient_threshold);
    ind_grad_high_high = find(train_dist_high(:,5) > gradient_threshold);
    
    train_bin1 = train_dist_low(ind_grad_low_low, :);  % Low dist + Low grad
    train_bin2 = train_dist_low(ind_grad_high_low, :);  % Low dist + High grad
    train_bin3 = train_dist_high(ind_grad_low_high, :); % High dist + Low grad
    train_bin4 = train_dist_high(ind_grad_high_high, :);% High dist + High grad
    
    % Check if all bins have sufficient samples (skip invalid splits)
    bin_sizes = [size(train_bin1,1), size(train_bin2,1), size(train_bin3,1), size(train_bin4,1)];
    
    if any(bin_sizes < min_bin_samples)
        if mod(total_attempts, 100) == 0
            fprintf('Attempt %d: Insufficient bin samples [%d, %d, %d, %d] - retrying...\n', ...
                total_attempts, bin_sizes(1), bin_sizes(2), bin_sizes(3), bin_sizes(4));
        end
        continue;  % Skip to next attempt if any bin is too small
    end
    
    % Valid split: increment counter and store metadata
    num_valid_iterations = num_valid_iterations + 1;
    
    % Build test matrix (zero out test areas to isolate train data)
    test_matrix = weight_matrix_40area;
    test_matrix(:, test_area_indices) = 0;  % Zero out test columns
    test_matrix(test_area_indices, :) = 0;  % Zero out test rows
    
    % Store train-test split metadata
    train_test_matrices{num_valid_iterations, 1} = train_matrix;   % Training matrix (91x91)
    train_test_matrices{num_valid_iterations, 2} = test_matrix;    % Testing matrix (91x91)
    train_test_matrices{num_valid_iterations, 3} = train_area_indices;  % Train area indices (20)
    train_test_matrices{num_valid_iterations, 4} = test_area_indices;   % Test area indices (20)
    
    % Train linear regression models for 4 bins (using global normalization)
    train_bins = {train_bin1, train_bin2, train_bin3, train_bin4};
    
    for bin_idx = 1:4
        bin_data = train_bins{bin_idx};
        
        % Normalize features (consistent with full data normalization)
        bin_data(:,8) = (bin_data(:,5) - nomal(bin_idx,1)) / nomal(bin_idx,2);  % Normalized gradient
        bin_data(:,9) = (bin_data(:,6) - nomal(bin_idx,3)) / nomal(bin_idx,4);  % Normalized distance
        
        % Train linear regression model (log(weight) ~ normalized gradient + distance)
        regression_models{num_valid_iterations, bin_idx} = fitlm(bin_data(:,8:9), bin_data(:,3));
        
        % Extract model coefficients (intercept + gradient + distance)
        intercept = table2array(regression_models{num_valid_iterations, bin_idx}.Coefficients(1,1));
        grad_coeff = table2array(regression_models{num_valid_iterations, bin_idx}.Coefficients(2,1));
        dist_coeff = table2array(regression_models{num_valid_iterations, bin_idx}.Coefficients(3,1));
        
        % Predict log weights for full bin data (all edges in the bin)
        predicted_log_weights{num_valid_iterations, bin_idx}(:,1) = intercept + ...
            full_bin_data{bin_idx}(:,8)*grad_coeff + full_bin_data{bin_idx}(:,9)*dist_coeff;
    end
    
    % Combine predictions from all 4 bins (full edge set)
    predicted_log_weights{num_valid_iterations, 5} = [predicted_log_weights{num_valid_iterations, 1}; ...
                                                      predicted_log_weights{num_valid_iterations, 2}; ...
                                                      predicted_log_weights{num_valid_iterations, 3}; ...
                                                      predicted_log_weights{num_valid_iterations, 4}];
    
    % Progress update every 100 valid iterations
    if mod(num_valid_iterations, 100) == 0
        fprintf('Completed %d/1000 valid iterations (total attempts: %d)\n', num_valid_iterations, total_attempts);
    end
end

% Final status update for train-test splits
if num_valid_iterations < 1000
    warning('Reached max attempts (%d) - only completed %d valid iterations', max_attempts, num_valid_iterations);
else
    fprintf('Successfully completed 1000 valid iterations (total attempts: %d)\n', total_attempts);
end

% Save intermediate results (for reproducibility)

save('matrain_train_test1000_40area_40_40.mat', 'train_test_matrices');
save('data4_40area_40_40.mat', 'full_bin_data');
save('data4_final_40area_40_40.mat', 'predicted_log_weights');

%% ===================== Part 3: ROC/AUC Calculation (Train/Test) =====================
% Reload saved results (optional - skip if running sequentially)

load('matrain_train_test1000_40area_40_40.mat');
load('data4_final_40area_40_40.mat');
load('MatF_SLN2021.mat');
weight_matrix_91x91 = MatF;

% Calculate ROC/AUC for TEST sets (1000 iterations)
test_roc_results = cell(1000, 2);
test_edge_counts = cell(1000, 1);

for iter_idx = 1:num_valid_iterations  % Use num_valid_iterations (handle incomplete runs)
    [test_roc_results{iter_idx,1}, test_roc_results{iter_idx,2}, test_edge_counts{iter_idx,1}] = ...
        cal_roc(predicted_log_weights{iter_idx,1}, train_test_matrices{iter_idx,2}, full_bin_data, weight_matrix_91x91, 20);
end

% Calculate ROC/AUC for TRAIN sets (1000 iterations)
fprintf('Processing TRAIN sets...\n');
train_roc_results = cell(1000, 2);
train_edge_counts = cell(1000, 1);

for iter_idx = 1:num_valid_iterations
    [train_roc_results{iter_idx,1}, train_roc_results{iter_idx,2}, train_edge_counts{iter_idx,1}] = ...
        cal_roc(predicted_log_weights{iter_idx,1}, train_test_matrices{iter_idx,1}, full_bin_data, weight_matrix_91x91, 20);
end

%% ===================== Part 4: ROC/AUC for 40-Area Full Prediction =====================
% Load full 40-area prediction results
load('data4_F_40area_40_40.mat');

% Calculate ROC/AUC for full 40-area prediction (all areas included)
[full_roc_sensitivity, full_roc_specificity, full_edge_details] = ...
    cal_roc(data4_F(:,10), weight_matrix_91x91, full_bin_data, weight_matrix_91x91, 40);

% Plot ROC curve for full 40-area prediction (publication-quality figure)
figure('Name', '40-Area ROC Curve', 'Position', [100, 100, 800, 600]);
plot(1 - full_roc_specificity, full_roc_sensitivity, 'r-', 'LineWidth', 1.5);
hold on;
line([0, 1], [0, 1], 'Color', 'gray', 'LineStyle', '--', 'LineWidth', 1);  % Random baseline
plot(1 - full_roc_specificity(5252,1), full_roc_sensitivity(5252,1), 'ro', 'MarkerSize', 6);
xlabel('1 - Specificity (False Positive Rate)', 'FontSize', 12);
ylabel('Sensitivity (True Positive Rate)', 'FontSize', 12);
title('ROC Curve: 40-Area Neuron Weight Prediction', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 10);
hold off;

% Calculate AUC (Area Under Curve) for full prediction
full_roc_specificity(:,2) = 1 - full_roc_specificity(:,1);
full_roc_diff = zeros(size(full_roc_specificity,1), 1);

for idx = 2:size(full_roc_specificity,1)
    full_roc_diff(idx) = full_roc_specificity(idx,2) - full_roc_specificity(idx-1,2);
end
full_roc_diff(1) = full_roc_specificity(1,2);  % First point difference (no prior point)

auc_full = trapz(full_roc_specificity(:,2), full_roc_sensitivity(:,1));


%% ===================== Part 5: Test Set Metrics (AUC/Sensitivity) =====================

% Process test set ROC data for AUC calculation
test_roc_processed = cell(num_valid_iterations, 1);
test_auc = zeros(num_valid_iterations, 1);
test_sensitivity = zeros(num_valid_iterations, 1);

for iter_idx = 1:num_valid_iterations
    % Convert specificity to 1-specificity (FPR)
    test_roc_processed{iter_idx} = test_roc_results{iter_idx,2};
    test_roc_processed{iter_idx}(:,2) = 1 - test_roc_processed{iter_idx}(:,1);
    
    % Calculate differences for AUC integration (trapezoidal rule)
    test_roc_processed{iter_idx}(:,3) = 0;
    if test_roc_processed{iter_idx}(1,2) > 0
        test_roc_processed{iter_idx}(1,3) = test_roc_processed{iter_idx}(1,2);
    end
    
    for idx = 2:8190
        test_roc_processed{iter_idx}(idx,3) = test_roc_processed{iter_idx}(idx,2) - test_roc_processed{iter_idx}(idx-1,2);
    end
    
    % Find index for 2274 edges (critical threshold) and calculate metrics
    edge_idx = find(test_edge_counts{iter_idx,1}(:,1)==2274, 1);
    if ~isempty(edge_idx)
        test_auc(iter_idx) = sum(test_roc_results{iter_idx,1}(:,1) .* test_roc_processed{iter_idx}(:,3));
        test_sensitivity(iter_idx) = test_roc_results{iter_idx,1}(edge_idx, 1);
    else
        test_auc(iter_idx) = NaN;
        test_sensitivity(iter_idx) = NaN;
        warning('Iteration %d: No 2274 edge index found', iter_idx);
    end
end

% Remove NaN values (if any)
test_auc = test_auc(~isnan(test_auc));
test_sensitivity = test_sensitivity(~isnan(test_sensitivity));

% Plot test set AUC distribution (publication-quality figure)
figure('Name', 'Test Set AUC Distribution', 'Position', [200, 200, 800, 600]);
hist(test_auc, 13);
xlabel('AUC', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Distribution of Test Set AUC (1000 Iterations)', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 10);
fprintf('Mean test set AUC: %.4f\n', mean(test_auc));

% Plot test set sensitivity distribution
figure('Name', 'Test Set Sensitivity Distribution', 'Position', [300, 300, 800, 600]);
hist(test_sensitivity, 13);
xlabel('Sensitivity (2274 Edges)', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Distribution of Test Set Sensitivity (1000 Iterations)', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 10);
fprintf('Mean test set sensitivity: %.4f\n', mean(test_sensitivity));

%% ===================== Part 6: Train Set Metrics (AUC/Sensitivity) =====================

% Process train set ROC data for AUC calculation
train_roc_processed = cell(num_valid_iterations, 1);
train_auc = zeros(num_valid_iterations, 1);
train_sensitivity = zeros(num_valid_iterations, 1);

for iter_idx = 1:num_valid_iterations
    % Convert specificity to 1-specificity (FPR)
    train_roc_processed{iter_idx} = train_roc_results{iter_idx,2};
    train_roc_processed{iter_idx}(:,2) = 1 - train_roc_processed{iter_idx}(:,1);
    
    % Calculate differences for AUC integration
    train_roc_processed{iter_idx}(:,3) = 0;
    if train_roc_processed{iter_idx}(1,2) > 0
        train_roc_processed{iter_idx}(1,3) = train_roc_processed{iter_idx}(1,2);
    end
    
    for idx = 2:8190
        train_roc_processed{iter_idx}(idx,3) = train_roc_processed{iter_idx}(idx,2) - train_roc_processed{iter_idx}(idx-1,2);
    end
    
    % Find index for 2274 edges and calculate metrics
    edge_idx = find(train_edge_counts{iter_idx,1}(:,1)==2274, 1);
    if ~isempty(edge_idx)
        train_auc(iter_idx) = sum(train_roc_results{iter_idx,1}(:,1) .* train_roc_processed{iter_idx}(:,3));
        train_sensitivity(iter_idx) = train_roc_results{iter_idx,1}(edge_idx, 1);
    else
        train_auc(iter_idx) = NaN;
        train_sensitivity(iter_idx) = NaN;
        warning('Iteration %d: No 2274 edge index found', iter_idx);
    end
end

% Remove NaN values (if any)
train_auc = train_auc(~isnan(train_auc));
train_sensitivity = train_sensitivity(~isnan(train_sensitivity));


%% ===================== Part 7: Random Baseline Analysis =====================
matrix_size = 91;
total_areas = 91;

% Extract actual non-zero edges (ground truth) for comparison
clear actual_edges;
actual_edge_count = 0;

for source_area = 1:total_areas
    for target_area = 1:total_areas
        if weight_matrix_91x91(source_area, target_area) ~= 0
            actual_edge_count = actual_edge_count + 1;
            actual_edges(actual_edge_count, 1) = source_area;
            actual_edges(actual_edge_count, 2) = target_area;
            actual_edges(actual_edge_count, 3) = weight_matrix_91x91(source_area, target_area);
        end
    end
end

% Find zero-degree columns (areas with no connections)
degree_col = sum(weight_matrix_91x91);
zero_degree_cols = find(degree_col == 0);

% Initialize random baseline storage
random_roc_results = cell(1000, 1);
random_edge_counts = cell(1000, 1);
random_auc = zeros(1000, 1);
random_sensitivity = zeros(1000, 1);

% 1000 iterations of random edge prediction (baseline comparison)
for rand_iter = 1:1000
    % Initialize random adjacency matrix (91x91)
    random_adj_matrix = zeros(matrix_size);
    
    % Randomly place edges (all non-self connections: 91×90=8190 edges)
    for edge_idx = 1:(matrix_size * (matrix_size - 1))
        % Find random non-self, non-duplicate edge
        while true
            rand_row = randi(matrix_size);
            rand_col = randi(matrix_size);
            if rand_row ~= rand_col && random_adj_matrix(rand_row, rand_col) == 0
                break;
            end
        end
        random_adj_matrix(rand_row, rand_col) = 1;
    end
    
    % Zero out zero-degree columns (match model constraints)
    filtered_rand_matrix = random_adj_matrix;
    for col_idx = 1:size(zero_degree_cols,1)
        filtered_rand_matrix(:, zero_degree_cols(col_idx)) = 0;
    end
    
    % Calculate TP/FP metrics (True Positive / Total Predicted)
    true_positive = 0;
    total_predicted = 0;
    actual_positive = size(actual_edges, 1);
    total_possible = (91*40) - 40;  % 40-area total possible edges
    actual_negative = total_possible - actual_positive;
    
    for source_area = 1:total_areas
        for target_area = 1:total_areas
            if filtered_rand_matrix(source_area, target_area) ~= 0
                total_predicted = total_predicted + 1;
                % Check if predicted edge is in ground truth
                if weight_matrix_91x91(source_area, target_area) ~= 0
                    true_positive = true_positive + 1;
                end
            end
        end
    end
    
    % Calculate sensitivity and specificity for random baseline
    if actual_positive > 0
        sensitivity = true_positive / actual_positive;
    else
        sensitivity = 0;
    end
    false_positive = total_predicted - true_positive;
    if actual_negative > 0
        specificity = (actual_negative - false_positive) / actual_negative;
    else
        specificity = 0;
    end
    
    % Store random ROC metrics
    random_roc_results{rand_iter,1}(edge_idx, 1) = sensitivity;
    random_roc_results{rand_iter,1}(edge_idx, 2) = specificity;
    random_edge_counts{rand_iter,1}(edge_idx, 1) = total_predicted;
end

% Process random baseline ROC data for AUC calculation
for rand_iter = 1:1000
    % Convert specificity to 1-specificity (FPR)
    random_roc_results{rand_iter,1}(:,3) = 1 - random_roc_results{rand_iter,1}(:,2);
    
    % Calculate differences for AUC integration
    random_roc_results{rand_iter,1}(:,4) = 0;
    if random_roc_results{rand_iter,1}(1,3) > 0
        random_roc_results{rand_iter,1}(1,4) = random_roc_results{rand_iter,1}(1,3);
    end
    
    for idx = 2:8190
        random_roc_results{rand_iter,1}(idx,4) = random_roc_results{rand_iter,1}(idx,3) - random_roc_results{rand_iter,1}(idx-1,3);
    end
    
    % Find index for 2274 edges and calculate metrics
    edge_idx = find(random_edge_counts{rand_iter,1}(:,1)==2274, 1);
    if ~isempty(edge_idx)
        random_auc(rand_iter) = sum(random_roc_results{rand_iter,1}(:,1) .* random_roc_results{rand_iter,1}(:,4));
        random_sensitivity(rand_iter) = random_roc_results{rand_iter,1}(edge_idx, 1);
    else
        random_auc(rand_iter) = NaN;
        random_sensitivity(rand_iter) = NaN;
    end
end

% Remove NaN values (if any)
random_auc = random_auc(~isnan(random_auc));
random_sensitivity = random_sensitivity(~isnan(random_sensitivity));

