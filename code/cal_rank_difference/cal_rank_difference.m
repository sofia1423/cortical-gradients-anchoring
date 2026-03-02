% =========================================================================
% @file    rank_difference_analysis.m
% @brief   Calculate empirical/predicted rank differences for spectral gradient data
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [Your Email]
% @note    Compare rank differences between empirical (467 samples) and predicted (91 regions) spectral gradients
% =========================================================================

%% ==================== Part 1: Empirical Rank Difference Calculation ====================
% Clear workspace and command window for clean execution
clearvars; clc; close all;

% Load empirical spectral gradient data (anesthesia close state)
load('Empirical_Spectral_Gradient_during_Anesthesia_and_Awake.mat');

% Flip sign of first eigenvector (consistent with original analysis)
empirical_gradient = U_ane_close;
empirical_gradient(:,1) = empirical_gradient(:,1) * (-1);

% Load region labels (467 samples × region indices)
load('Label_467.mat');

% Initialize empirical data matrix (467 samples)
empirical_data = zeros(467, 7);
empirical_data(:,1) = 1:467;                          % Sample indices (1-467)
empirical_data(:,2:3) = empirical_gradient(:,1:2);     % First two spectral gradient components

% Calculate ranking for each gradient component (descending order)
for comp_idx = 1:2
    % Extract current gradient component
    gradient_component = empirical_data(:, comp_idx + 1);
    
    % Sort component in descending order and get sorted indices
    [sorted_values, sorted_indices] = sort(gradient_component, 'descend');
    
    % Initialize ranking vector
    component_ranking = zeros(size(gradient_component));
    
    % Assign ranks (1 = highest value)
    for sample_idx = 1:length(gradient_component)
        component_ranking(sorted_indices(sample_idx)) = sample_idx;
    end
    
    % Store ranking in empirical data matrix (columns 4,5)
    empirical_data(:, 3 + comp_idx) = component_ranking;
end

% Calculate rank difference (component 2 rank - component 1 rank)
empirical_data(:,6) = empirical_data(:,5) - empirical_data(:,4);

% Normalize rank difference (z-score)
empirical_data(:,6) = zscore(empirical_data(:,6));

% Add region labels to empirical data
empirical_data(:,7) = label_467;

% Group samples by region (91 total regions, filter non-empty regions)
region_groups = cell(91, 1);
valid_region_indices = [];
valid_region_count = 0;

for region_idx = 1:91
    % Find samples belonging to current region
    region_sample_indices = find(empirical_data(:,7) == region_idx);
    
    if ~isempty(region_sample_indices)
        valid_region_count = valid_region_count + 1;
        valid_region_indices(valid_region_count, 1) = region_idx;
        
        % Store region sample indices and their z-scored rank differences
        region_groups{region_idx}(:,1) = region_sample_indices;
        region_groups{region_idx}(:,2) = empirical_data(region_sample_indices, 6);
    end
end

% Filter to only valid (non-empty) regions
region_groups = region_groups(valid_region_indices);

% Calculate mean rank difference per region (empirical)
region_rank_diff_empirical = cell(valid_region_count, 1);
empirical_rank_difference = zeros(91, 1);

for region_idx = 1:valid_region_count
    % Extract z-scored rank differences for current region
    region_rank_diff_empirical{region_idx} = region_groups{region_idx}(:,2);
    
    % Calculate mean rank difference for the region
    empirical_mean = mean(region_rank_diff_empirical{region_idx});
    empirical_rank_difference(valid_region_indices(region_idx)) = empirical_mean;
end

%% ==================== Part 2: Predicted Rank Difference Calculation ====================
% Load predicted spectral gradient data (anesthesia + awake states)
load('Predicted_Spectral_Gradient_during_Anesthesia_and_Awake.mat');

% Initialize predicted data matrix (91 regions)
predicted_data = zeros(91, 6);
predicted_data(:,1) = 1:91;                          % Region indices (1-91)
predicted_data(:,2:3) = spectral_gradient(:,1:2);    % First two predicted gradient components

% Calculate ranking for each predicted component (ascending order - match original)
for comp_idx = 1:2
    % Extract current predicted component
    pred_component = predicted_data(:, comp_idx + 1);
    
    % Sort component in ascending order and get sorted indices
    [sorted_pred_values, sorted_pred_indices] = sort(pred_component, 'ascend');
    
    % Initialize ranking vector
    pred_component_ranking = zeros(size(pred_component));
    
    % Assign ranks (1 = lowest value)
    for region_idx = 1:length(pred_component)
        pred_component_ranking(sorted_pred_indices(region_idx)) = region_idx;
    end
    
    % Store ranking in predicted data matrix (columns 4,5)
    predicted_data(:, 3 + comp_idx) = pred_component_ranking;
end

% Calculate rank difference (component 2 rank - component 1 rank)
predicted_data(:,6) = predicted_data(:,5) - predicted_data(:,4);

% Normalize rank difference (z-score)
predicted_rank_difference = zscore(predicted_data(:,6));
