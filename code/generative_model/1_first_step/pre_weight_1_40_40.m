% =========================================================================
% @file    neuron_connectivity_binning_fitting.m
% @brief   Neuron Connectivity Analysis: 2x2 Binning & Linear Regression Fitting
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key workflow: 
%          1. Preprocess 91-area weight matrix to 40x40 valid area matrix
%          2. Build feature matrix with gradient/distance/weight features
%          3. Sweep distance (0.1-70) & gradient (0.01-1) thresholds for 2x2 binning
%          4. Fit linear models to each bin & filter results for negative coefficients
% =========================================================================

%% ===================== Part 1: Data Loading & Preprocessing =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;
rng('shuffle');  % Set random seed (for any stochastic operations)

% Load core datasets (critical dependencies)
fprintf('Loading core datasets...\n');
load('MatF_SLN2021.mat');          % 91x91 neuron weight matrix (ground truth)
load('lp19_gd_fitnew.mat');        % Gradient fit results (gdfit10: 91x1 gradient values)
load('distance_matrix.mat');       % 91x91 pairwise distance matrix (dist_m)

% Preprocess weight matrix: remove zero columns (filter valid brain areas)
zero_column_indices = find(all(MatF == 0, 1));  % Columns with all zero values (invalid areas)
full_area_indices = 1:91;                       % Initial 91-area index range
valid_area_indices_91 = full_area_indices;
valid_area_indices_91(zero_column_indices) = [];% Remove zero columns (valid 40 areas)

% Transform to 40x40 weight matrix (using custom trans_matrix_40area function)
weight_matrix_40x40 = trans_matrix_40area(MatF, 40, 40);

%% ===================== Part 2: Build Neuron Feature Matrix =====================
% Initialize feature matrix to store core connectivity features
% Columns: 1=row_idx, 2=col_idx, 3=log_weight, 4=raw_weight, 5=abs_grad_ratio, 6=distance, 10=mapped_row, 11=mapped_col
clear neuron_feature_matrix;
feature_row_count = 0;
matrix_size = size(weight_matrix_40x40);

% Iterate through non-zero entries in 40x40 weight matrix (vectorized for speed)

non_zero_indices = find(weight_matrix_40x40 ~= 0);  % Vectorized non-zero lookup
[weight_row_indices, weight_col_indices] = ind2sub(matrix_size, non_zero_indices);

for idx = 1:length(non_zero_indices)
    weight_row = weight_row_indices(idx);
    weight_col = weight_col_indices(idx);
    feature_row_count = feature_row_count + 1;
    
    % 1. Core index mapping (40-area to original 91-area indices)
    neuron_feature_matrix(feature_row_count, 1)  = weight_row;                  % Row index (40x40 matrix)
    neuron_feature_matrix(feature_row_count, 2)  = weight_col;                  % Column index (40x40 matrix)
    neuron_feature_matrix(feature_row_count, 10) = valid_area_indices_91(weight_row); % Mapped 91-area row index
    neuron_feature_matrix(feature_row_count, 11) = valid_area_indices_91(weight_col); % Mapped 91-area column index
    
    % 2. Weight features (raw + log-transformed)
    raw_weight = weight_matrix_40x40(weight_row, weight_col);
    neuron_feature_matrix(feature_row_count, 3)  = log(raw_weight);             % Log-transformed weight (for regression)
    neuron_feature_matrix(feature_row_count, 4)  = raw_weight;                  % Raw weight value
    
    % 3. Gradient & distance features (using mapped 91-area indices)
    grad_ratio = log(gdfit10(valid_area_indices_91(weight_row), 1) / ...
                     gdfit10(valid_area_indices_91(weight_col), 1));
    neuron_feature_matrix(feature_row_count, 5)  = abs(grad_ratio);             % Absolute log gradient ratio
    neuron_feature_matrix(feature_row_count, 6)  = dist_m(valid_area_indices_91(weight_row), ...
                                                          valid_area_indices_91(weight_col)); % Pairwise distance
end

%% ===================== Part 3: 2x2 Binning & Linear Fitting Pipeline =====================
% Initialize output matrix for binning/fitting results
% Column layout:
% 1-4: Bin sample counts (Bin1/Bin2/Bin3/Bin4)
% 5-8: Sqrt(R²) for each bin's regression model
% 9-10: Global correlation coeff + p-value (fitted vs original log weights)
% 11-12: Distance/gradient thresholds used for binning
% 13+: Regression coefficients (grad/dist for each bin)
clear fitting_results; 
valid_bin_count = 0;

% Define threshold sweep parameters (vectorized for speed)
distance_thresholds = 0.1:0.1:70;    % Distance sweep (0.1 to 70, step 0.1: 700 values)
gradient_thresholds = 0.01:0.01:1;  % Gradient sweep (0.01 to 1, step 0.01: 100 values)

% Precompute feature columns for faster access
grad_feature = neuron_feature_matrix(:,5);   % Absolute gradient ratio (column 5)
dist_feature = neuron_feature_matrix(:,6);   % Pairwise distance (column 6)
log_weight_feature = neuron_feature_matrix(:,3); % Log-transformed weight (column 3)

% Main threshold sweep loop (distance × gradient)
progress_counter = 0;

for dist_idx = 1:length(distance_thresholds)
    distance_threshold = distance_thresholds(dist_idx);
    
    % Split feature matrix by distance threshold (vectorized logical indexing)
    is_low_distance = dist_feature <= distance_threshold;
    feature_low_distance = neuron_feature_matrix(is_low_distance, :);
    feature_high_distance = neuron_feature_matrix(~is_low_distance, :);
    
    % Extract gradient features for low/high distance subsets (speed optimization)
    grad_low_dist = feature_low_distance(:,5);
    grad_high_dist = feature_high_distance(:,5);

    for grad_idx = 1:length(gradient_thresholds)
        progress_counter = progress_counter + 1;
        gradient_threshold = gradient_thresholds(grad_idx);
        
        % Split low/high distance data by gradient threshold (vectorized)
        is_low_grad_low_dist = grad_low_dist <= gradient_threshold;
        is_high_grad_low_dist = grad_low_dist > gradient_threshold;
        is_low_grad_high_dist = grad_high_dist <= gradient_threshold;
        is_high_grad_high_dist = grad_high_dist > gradient_threshold;
        
        % 2x2 binning (4 bins: distance × gradient combinations)
        bin1 = feature_low_distance(is_low_grad_low_dist, :);   % Bin1: Low dist + Low grad
        bin2 = feature_low_distance(is_high_grad_low_dist, :);  % Bin2: Low dist + High grad
        bin3 = feature_high_distance(is_low_grad_high_dist, :); % Bin3: High dist + Low grad
        bin4 = feature_high_distance(is_high_grad_high_dist, :);% Bin4: High dist + High grad

        % Only process if all bins have ≥2 samples (statistical validity for regression)
        bin_sizes = [size(bin1,1), size(bin2,1), size(bin3,1), size(bin4,1)];
        if all(bin_sizes > 1)
            valid_bin_count = valid_bin_count + 1;
            
            % Store bin metadata (thresholds + sample counts)
            fitting_results(valid_bin_count, 11) = distance_threshold;  % Distance split threshold
            fitting_results(valid_bin_count, 12) = gradient_threshold;  % Gradient split threshold
            fitting_results(valid_bin_count, 1:4) = bin_sizes;           % Sample count per bin

            % Initialize arrays for global correlation calculation
            fitted_log_weights = [];  % Fitted log-weight values from regression
            original_log_weights = [];% Original log-weight values from feature matrix

            % Process each bin (normalization + linear regression)
            bin_list = {bin1, bin2, bin3, bin4};
            coeff_col_offset = 12;  % Starting column for regression coefficients
            
            for bin_index = 1:4
                bin_data = bin_list{bin_index};
                
                % Min-Max normalization (scale features to [0,1] for regression stability)
                grad_min = min(bin_data(:,5));
                grad_max = max(bin_data(:,5));
                dist_min = min(bin_data(:,6));
                dist_max = max(bin_data(:,6));
                
                % Avoid division by zero (edge case: all values identical)
                bin_data(:,8) = (bin_data(:,5) - grad_min) / (grad_max - grad_min + eps); % Normalized gradient
                bin_data(:,9) = (bin_data(:,6) - dist_min) / (dist_max - dist_min + eps); % Normalized distance

                % Linear regression: log(weight) ~ normalized(gradient + distance)
                regression_model = fitlm(bin_data(:,8:9), bin_data(:,3));

                % Store regression performance (sqrt of R² for interpretability)
                fitting_results(valid_bin_count, bin_index + 4) = sqrt(regression_model.Rsquared.Ordinary);
                
                % Extract & store regression coefficients (gradient + distance)
                coeff_table = regression_model.Coefficients;
                grad_coeff = table2array(coeff_table(2,1));  % Gradient feature coefficient
                dist_coeff = table2array(coeff_table(3,1));  % Distance feature coefficient
                
                fitting_results(valid_bin_count, coeff_col_offset + (bin_index-1)*2 + 1) = grad_coeff;
                fitting_results(valid_bin_count, coeff_col_offset + (bin_index-1)*2 + 2) = dist_coeff;

                % Accumulate fitted/original values for global correlation
                fitted_log_weights = [fitted_log_weights; regression_model.Fitted];
                original_log_weights = [original_log_weights; bin_data(:,3)];
            end

            % Calculate global correlation between fitted & original log weights
            if ~isempty(fitted_log_weights) && ~isempty(original_log_weights)
                [corr_coeff, p_value] = corr_lp(fitted_log_weights, original_log_weights);
                fitting_results(valid_bin_count, 9) = corr_coeff;  % Correlation coefficient
                fitting_results(valid_bin_count, 10) = p_value;   % P-value of correlation
            else
                fitting_results(valid_bin_count, 9:10) = NaN;  % Handle empty data case
            end
        end
        
        % Progress update every 1000 iterations
        if mod(progress_counter, 1000) == 0
            fprintf('Progress: %d/%d iterations completed\n', progress_counter, length(distance_thresholds)*length(gradient_thresholds));
        end
    end
end

%% ===================== Part 4: Filter Results (Negative Coefficients Only) =====================
% Define column indices for regression coefficients (4 bins × 2 features)
coeff_columns = [13,14, 16,17, 19,20, 22,23];  % Grad/Dist for Bin1-Bin4

% Filter rows where ALL 8 coefficients are negative (strict negative filter)
is_negative_coeff = true(size(fitting_results,1),1);
for col_idx = coeff_columns
    is_negative_coeff = is_negative_coeff & (fitting_results(:,col_idx) < 0);
end

% Extract final results (only negative coefficient threshold combinations)
negative_coeff_indices = find(is_negative_coeff);
final_results_negative_coeff = fitting_results(negative_coeff_indices, :);
