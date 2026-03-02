% =========================================================================
% @file    receptor_extraction_svd.m
% @brief   Extract receptor eigenvector from transmitter data using SVD
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [Your Email]
% @note    Process 91-region transmitter data to get 1st SVD eigenvector (receptor)
% =========================================================================

%% ==================== Data Loading & Preprocessing ====================
% Clear workspace and command window for clean execution
clearvars; clc; close all;

% Load transmitter data (91 regions × 14 features)
load('transmitter_1_14.mat');
transmitter_data = transmitter;  % Rename loaded variable for clarity

% Filter non-zero transmitter regions (remove zero rows)
valid_region_indices = [];
filtered_transmitter_data = [];
valid_region_count = 0;

for region_idx = 1:91
    if transmitter_data(region_idx, 1) > 0  % Keep regions with non-zero first feature
        valid_region_count = valid_region_count + 1;
        valid_region_indices(valid_region_count, 1) = region_idx;          % Store valid region indices
        filtered_transmitter_data(valid_region_count, :) = transmitter_data(region_idx, :);  % Filtered data
    end
end

fprintf('Filtered transmitter data: %d valid regions (from 91 total)\n', valid_region_count);

%% ==================== Data Normalization (Mean Centering) ====================
num_features = 14;  % 14 transmitter features (hardcoded as per original code)
feature_means = zeros(1, num_features);  % Mean of each feature across valid regions
mean_centered_data = zeros(size(filtered_transmitter_data));  % Mean-centered data

% Calculate feature means and center data (subtract mean)
for feature_idx = 1:num_features
    feature_means(1, feature_idx) = mean(filtered_transmitter_data(:, feature_idx));
    mean_centered_data(:, feature_idx) = filtered_transmitter_data(:, feature_idx) - feature_means(1, feature_idx);
end

%% ==================== L2 Normalization (Unit Norm) ====================
% Calculate total sum of squares for L2 normalization
total_squared_sum = 0;
for row_idx = 1:size(mean_centered_data, 1)
    for col_idx = 1:size(mean_centered_data, 2)
        total_squared_sum = total_squared_sum + mean_centered_data(row_idx, col_idx)^2;
    end
end

% Compute L2 norm and normalize data to unit norm
l2_norm = sqrt(total_squared_sum);
normalized_data = mean_centered_data / l2_norm;  % Final normalized matrix (A) for SVD

fprintf('Data normalized to unit L2 norm (norm value: %.4f)\n', l2_norm);

%% ==================== Singular Value Decomposition (SVD) ====================
% Perform SVD: normalized_data = U * S * V'
% - U: Left singular vectors (valid regions × valid regions)
% - S: Singular values (diagonal matrix)
% - V: Right singular vectors (features × features)
[U, S, V] = svd(normalized_data);

% Extract 1st eigenvector from U (receptor final)
receptor_eigenvector = U(:, 1);

fprintf('Receptor eigenvector extracted (length: %d)\n', length(receptor_eigenvector));