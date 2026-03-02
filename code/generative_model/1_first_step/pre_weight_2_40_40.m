% =========================================================================
% @file    neuron_edge_weight_prediction.m
% @brief   Neuron Edge Weight Prediction (8190 Edges): 4-Bin Linear Model Pipeline
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key workflow:
%          1. Preprocess 91-area weight matrix to 40x40 valid area matrix
%          2. Train 4-bin linear models (distance=33.5, gradient=0.19 thresholds)
%          3. Predict weights for all 8190 non-self edges in 91x91 matrix
%          4. Validate predictions against training data (non-zero weights)
% =========================================================================

%% ===================== Part 1: Environment Setup & Data Loading =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;
rng('shuffle');  % Set random seed (for consistent regression results)

% Load core datasets (critical dependencies)
fprintf('Loading core datasets...\n');
load('MatF_SLN2021.mat');               % 91x91 neuron weight matrix (ground truth)
load('lp19_gd_fitnew.mat');             % Gradient fit results (gdfit10: 91x1 gradient values)
load('distance_matrix.mat');            % 91x91 pairwise distance matrix (dist_m)

% Preprocess weight matrix: remove zero columns (filter valid brain areas)
zero_column_indices = find(all(MatF == 0, 1));  % Columns with all zero values (invalid areas)
full_area_indices = 1:91;                       % Initial 91-area index range
valid_area_indices_91 = full_area_indices;
valid_area_indices_91(zero_column_indices) = [];% Retain only valid areas (40 total)

% Transform to 40x40 weight matrix (using custom transformation function)
weight_matrix_40x40 = trans_matrix_40area(MatF, 40, 40);

%% ===================== Part 2: Build Training Feature Matrix (Non-Zero Weights) =====================
% Initialize feature matrix for training data (non-zero weight entries)
% Columns: 1=row_idx(40x40), 2=col_idx(40x40), 3=log_weight, 4=raw_weight, 5=abs_grad_ratio, 6=distance, 10=mapped_row, 11=mapped_col
clear training_feature_matrix;
feature_row_count = 0;
matrix_size = size(weight_matrix_40x40);

% Extract non-zero entries (vectorized for speed)

non_zero_indices = find(weight_matrix_40x40 ~= 0);
[weight_row_indices, weight_col_indices] = ind2sub(matrix_size, non_zero_indices);

for idx = 1:length(non_zero_indices)
    weight_row = weight_row_indices(idx);
    weight_col = weight_col_indices(idx);
    feature_row_count = feature_row_count + 1;
    
    % 1. Core index mapping (40x40 to original 91-area indices)
    training_feature_matrix(feature_row_count, 1)  = weight_row;                  % Row index (40x40 matrix)
    training_feature_matrix(feature_row_count, 2)  = weight_col;                  % Column index (40x40 matrix)
    training_feature_matrix(feature_row_count, 10) = valid_area_indices_91(weight_row); % Mapped 91-area row index
    training_feature_matrix(feature_row_count, 11) = valid_area_indices_91(weight_col); % Mapped 91-area column index
    
    % 2. Weight features (raw + log-transformed for regression)
    raw_weight = weight_matrix_40x40(weight_row, weight_col);
    training_feature_matrix(feature_row_count, 3)  = log(raw_weight);             % Log-transformed weight
    training_feature_matrix(feature_row_count, 4)  = raw_weight;                  % Raw weight value
    
    % 3. Gradient & distance features (using mapped 91-area indices)
    grad_ratio = log(gdfit10(valid_area_indices_91(weight_row), 1) / ...
                     gdfit10(valid_area_indices_91(weight_col), 1));
    training_feature_matrix(feature_row_count, 5)  = abs(grad_ratio);             % Absolute log gradient ratio
    training_feature_matrix(feature_row_count, 6)  = dist_m(valid_area_indices_91(weight_row), ...
                                                          valid_area_indices_91(weight_col)); % Pairwise distance
end


%% ===================== Part 3: Define 4 Bins (Thresholds for 40x40 Version) =====================
% Critical thresholds (2026.2.7 update: optimized for 40x40 weight matrix)
distance_threshold = 33.5;  % Distance cutoff (mm)
gradient_threshold = 0.19;  % Gradient ratio cutoff

% Assign training data to 4 bins (logical indexing for speed)
bin_indices{1} = find((training_feature_matrix(:,6) <= distance_threshold) & (training_feature_matrix(:,5) <= gradient_threshold));
bin_indices{2} = find((training_feature_matrix(:,6) > distance_threshold) & (training_feature_matrix(:,5) <= gradient_threshold));
bin_indices{3} = find((training_feature_matrix(:,6) <= distance_threshold) & (training_feature_matrix(:,5) > gradient_threshold));
bin_indices{4} = find((training_feature_matrix(:,6) > distance_threshold) & (training_feature_matrix(:,5) > gradient_threshold));

% Print bin sizes for verification
bin_sizes = [length(bin_indices{1}), length(bin_indices{2}), length(bin_indices{3}), length(bin_indices{4})];

%% ===================== Part 4: Train 4-Bin Linear Models (Normalization + Fitting) =====================
% Initialize storage for bin data, models, and normalization parameters
clear bin_data_list;
clear regression_models;
normalization_params = zeros(4,4);  % [grad_min, grad_range, dist_min, dist_range] per bin
fitted_log_weights = [];
original_log_weights = [];

for bin_idx = 1:4
    % Extract data for current bin
    bin_data_list{bin_idx} = training_feature_matrix(bin_indices{bin_idx}, :);
    
    % Min-Max normalization for gradient (column 5) - add eps to avoid division by zero
    grad_min = min(bin_data_list{bin_idx}(:,5));
    grad_max = max(bin_data_list{bin_idx}(:,5));
    grad_range = grad_max - grad_min + eps;
    bin_data_list{bin_idx}(:,8) = (bin_data_list{bin_idx}(:,5) - grad_min) / grad_range;
    normalization_params(bin_idx, 1:2) = [grad_min, grad_range];  % Store gradient normalization params
    
    % Min-Max normalization for distance (column 6)
    dist_min = min(bin_data_list{bin_idx}(:,6));
    dist_max = max(bin_data_list{bin_idx}(:,6));
    dist_range = dist_max - dist_min + eps;
    bin_data_list{bin_idx}(:,9) = (bin_data_list{bin_idx}(:,6) - dist_min) / dist_range;
    normalization_params(bin_idx, 3:4) = [dist_min, dist_range];  % Store distance normalization params
    
    % Train linear regression model: log(weight) ~ normalized(gradient + distance)
    regression_models{bin_idx} = fitlm(bin_data_list{bin_idx}(:,8:9), bin_data_list{bin_idx}(:,3));
    
    % Accumulate fitted/original values for training correlation
    fitted_log_weights = [fitted_log_weights; regression_models{bin_idx}.Fitted];
    original_log_weights = [original_log_weights; bin_data_list{bin_idx}(:,3)];
end

% Calculate training correlation (validation on non-zero weights)
[train_corr, train_p_value] = corr_lp(original_log_weights, fitted_log_weights);

% Save normalization parameters (for downstream use)
save('Normal40area_40_40.mat', 'normalization_params');

%% ===================== Part 5: Predict Weights for All 8190 Edges (Full 91x91 Matrix) =====================
% Build full edge matrix (all non-self edges: 91x90 = 8190 edges)

clear full_edge_matrix;
edge_count = 0;

for source_area = 1:91
    for target_area = 1:91
        if source_area ~= target_area  % Exclude self-connections (i==k)
            edge_count = edge_count + 1;
            full_edge_matrix(edge_count, 1) = source_area;          % Source area index (91x91)
            full_edge_matrix(edge_count, 2) = target_area;          % Target area index (91x91)
            full_edge_matrix(edge_count, 3) = abs(log(gdfit10(source_area,1)/gdfit10(target_area,1)));  % Gradient ratio
            full_edge_matrix(edge_count, 4) = dist_m(source_area, target_area);  % Pairwise distance
        end
    end
end

% Assign full edges to 4 bins (same thresholds as training)
full_bin_indices{1} = find((full_edge_matrix(:,4) <= distance_threshold) & (full_edge_matrix(:,3) <= gradient_threshold));
full_bin_indices{2} = find((full_edge_matrix(:,4) > distance_threshold) & (full_edge_matrix(:,3) <= gradient_threshold));
full_bin_indices{3} = find((full_edge_matrix(:,4) <= distance_threshold) & (full_edge_matrix(:,3) > gradient_threshold));
full_bin_indices{4} = find((full_edge_matrix(:,4) > distance_threshold) & (full_edge_matrix(:,3) > gradient_threshold));

% Predict log weights for each bin using trained models
clear full_bin_data;

for bin_idx = 1:4
    full_bin_data{bin_idx} = full_edge_matrix(full_bin_indices{bin_idx}, :);
    
    % Apply training-set normalization (critical for consistent predictions)
    full_bin_data{bin_idx}(:,8) = (full_bin_data{bin_idx}(:,3) - normalization_params(bin_idx,1)) / normalization_params(bin_idx,2);
    full_bin_data{bin_idx}(:,9) = (full_bin_data{bin_idx}(:,4) - normalization_params(bin_idx,3)) / normalization_params(bin_idx,4);
    
    % Extract regression coefficients (intercept + gradient + distance)
    coeff_table = regression_models{bin_idx}.Coefficients;
    intercept = table2array(coeff_table(1,1));
    grad_coeff = table2array(coeff_table(2,1));
    dist_coeff = table2array(coeff_table(3,1));
    
    % Predict log weights for current bin
    full_bin_data{bin_idx}(:,10) = intercept + full_bin_data{bin_idx}(:,8)*grad_coeff + full_bin_data{bin_idx}(:,9)*dist_coeff;
end

% Combine all predicted edges into a single matrix
predicted_full_edges = [full_bin_data{1}; full_bin_data{2}; full_bin_data{3}; full_bin_data{4}];

%% ===================== Part 6: Validate Predictions Against Training Data =====================

validation_edge_count = size(training_feature_matrix, 1);
prediction_indices = zeros(validation_edge_count, 1);

% Match predicted edges to training edges (by source/target area)
% Vectorized lookup for speed (replaces slow loop for large datasets)
for idx = 1:validation_edge_count
    source_area = training_feature_matrix(idx, 10);
    target_area = training_feature_matrix(idx, 11);
    match_idx = find((predicted_full_edges(:,1)==source_area) & (predicted_full_edges(:,2)==target_area), 1);
    if ~isempty(match_idx)
        prediction_indices(idx,1) = match_idx;
    else
        prediction_indices(idx,1) = NaN;  % Handle missing edges (edge case)
    end
end

% Remove NaN entries (if any)
valid_pred_indices = ~isnan(prediction_indices);
predicted_log_weights = predicted_full_edges(prediction_indices(valid_pred_indices), 10);
actual_log_weights = training_feature_matrix(valid_pred_indices, 3);

% Calculate validation correlation
if ~isempty(predicted_log_weights) && ~isempty(actual_log_weights)
    [val_corr, val_p_value] = corr_lp(actual_log_weights, predicted_log_weights);
    fprintf('Validation correlation (training edges): %.4f (p-value: %.4e)\n', val_corr, val_p_value);
else
    warning('No valid prediction matches found for validation');
end

%% ===================== Part 7: Build Full 91x91 Predicted Weight Matrix =====================
predicted_weight_matrix = zeros(91,91);  % Initialize empty 91x91 matrix

% Assign predicted log weights to 91x91 matrix
for idx = 1:size(predicted_full_edges, 1)
    source_area = predicted_full_edges(idx, 1);
    target_area = predicted_full_edges(idx, 2);
    predicted_weight_matrix(source_area, target_area) = predicted_full_edges(idx, 10);
end

% Save final predicted matrix (critical output for downstream analysis)
save('matrix_pre_8190_40area.mat', 'predicted_weight_matrix');
