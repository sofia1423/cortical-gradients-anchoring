% =========================================================================
% @file    neuron_gradient_binning_91x40.m
% @brief   Neuron Gradient Binning & Distribution Comparison (91x40 Area Matrix)
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key workflow:
%          1. Compute 91x91 gradient ratio matrix (log(gdfit10(k)/gdfit10(i)))
%          2. Transform to 91x40 matrix using custom trans_matrix_40area
%          3. Bin gradient values (variable bin count: 2*(ttt+2))
%          4. Compare bin counts between Data_B and Data_A
%          5. Calculate Pdif (difference between symmetric bins)
% =========================================================================

%% ===================== Part 1: Environment Setup & Initialization =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;
rng('shuffle');  % Set random seed for consistent results

% Define core parameters (13 iterations, bin count = 2*(ttt+2))
total_iterations = 13;
data_nd_618 = cell(total_iterations, 2);  % Store P_final (1) and thr_Bin (2)
Pdif = cell(total_iterations, 1);         % Store Pdif results (symmetric bin differences)

%% ===================== Part 2: Main Binning Loop (13 Iterations) =====================
for iteration_idx = 1:total_iterations
    fprintf('\n=== Processing iteration %d/%d ===\n', iteration_idx, total_iterations);
    
    % Calculate bin count for current iteration (2*(ttt+2))
    bin_count = 2 * (iteration_idx + 2);
    fprintf('Bin count for iteration %d: %d\n', iteration_idx, bin_count);
    
    %% Step 2.1: Load Core Data & Compute 91x91 Gradient Ratio Matrix
    load('lp19_gd_fitnew.mat');  % Gradient fit results (gdfit10: 91x1)
    load('MatF_SLN2021.mat');    % Neuron weight matrix (for 40-area transformation)
    
    % Initialize 91x91 gradient ratio matrix (log(gdfit10(k,1)/gdfit10(i,1)))
    gradient_ratio_matrix = zeros(91, 91);
    for source_area = 1:91
        for target_area = 1:91
            gradient_ratio_matrix(source_area, target_area) = log(gdfit10(target_area, 1) / gdfit10(source_area, 1));
        end
    end
  
    
    %% Step 2.2: Transform to 91x40 Matrix & Extract Feature Arrays
    % Transform gradient matrix to 91x40 (valid area matrix)
    gradient_matrix_91x40 = trans_matrix_40area(gradient_ratio_matrix, 91, 40);
    
    % Extract feature arrays (custom functions: data_array_40area/data_array2_40area)
    Data_B = data_array_40area(gradient_matrix_91x40);
    Data_A = data_array2_40area(gradient_matrix_91x40);
   
    %% Step 2.3: Normalize Data_B to Full Gradient Range
    % Get global min/max of original gradient matrix (for normalization)
    grad_min_global = min(gradient_ratio_matrix, [], 'all');
    grad_max_global = max(gradient_ratio_matrix, [], 'all');
    
    % Min-Max normalize Data_B column 3 to [grad_min_global, grad_max_global]
    Data_B_col3_min = min(Data_B(:,3));
    Data_B_col3_max = max(Data_B(:,3));
    Data_B(:,4) = (grad_max_global - grad_min_global) * (Data_B(:,3) - Data_B_col3_min) / ...
                  (Data_B_col3_max - Data_B_col3_min + eps) + grad_min_global;  % +eps to avoid division by zero
    
    %% Step 2.4: Rebuild 91x40 Matrix from Normalized Data_B
    gradient_matrix_91x40_norm = zeros(91, 40);  % Initialize normalized 91x40 matrix
    for row_idx = 1:size(Data_B,1)
        source_idx = Data_B(row_idx, 1);
        target_idx = Data_B(row_idx, 2);
        gradient_matrix_91x40_norm(source_idx, target_idx) = Data_B(row_idx, 4);
    end
    
    %% Step 2.5: Define Bin Thresholds & Count Samples per Bin
    % Calculate bin interval (evenly spaced between min/max of Data_B column 4)
    Data_B_col4_min = min(Data_B(:,4));
    Data_B_col4_max = max(Data_B(:,4));
    bin_interval = (Data_B_col4_max - Data_B_col4_min) / bin_count;
    
    % Initialize bin thresholds and count matrices
    bin_thresholds = zeros(bin_count, 2);  % [lower, upper] per bin
    bin_counts = zeros(2, bin_count);     % Row1=Data_B, Row2=Data_A
    
    % Calculate bin thresholds and count samples in each bin
    for bin_idx = 1:bin_count
        % Define bin boundaries (lower = min + (i-1)*interval, upper = min + i*interval)
        bin_thresholds(bin_idx, 1) = Data_B_col4_min + (bin_idx - 1) * bin_interval;
        bin_thresholds(bin_idx, 2) = Data_B_col4_min + bin_idx * bin_interval;
        
        % Count samples in current bin (Data_B column 4)
        bin_counts(1, bin_idx) = sum(Data_B(:,4) > bin_thresholds(bin_idx,1) & Data_B(:,4) <= bin_thresholds(bin_idx,2));
        
        % Count samples in current bin (Data_A column 3)
        bin_counts(2, bin_idx) = sum(Data_A(:,3) > bin_thresholds(bin_idx,1) & Data_A(:,3) <= bin_thresholds(bin_idx,2));
    end
    
    %% Step 2.6: Adjust Bin Counts (Ensure Total Matches Data Size)
    % Adjust Data_B counts if sum < total samples
    total_Data_B = size(Data_B, 1);
    sum_Data_B = sum(bin_counts(1,:));
    if sum_Data_B < total_Data_B
        bin_counts(1,1) = bin_counts(1,1) + (total_Data_B - sum_Data_B);  % Add to first bin
        fprintf('Adjusted Data_B bin 1 count: +%d (total mismatch)\n', total_Data_B - sum_Data_B);
    end
    
    % Adjust Data_A counts if sum < total samples
    total_Data_A = size(Data_A, 1);
    sum_Data_A = sum(bin_counts(2,:));
    if sum_Data_A < total_Data_A
        bin_counts(2,1) = bin_counts(2,1) + (total_Data_A - sum_Data_A);  % Add to first bin
        fprintf('Adjusted Data_A bin 1 count: +%d (total mismatch)\n', total_Data_A - sum_Data_A);
    end
    
    
    %% Step 2.7: Calculate P_final (Data_A count / Data_B count)
    % Avoid division by zero (replace 0 with eps in Data_B counts)
    bin_counts_Data_B_safe = bin_counts(1,:);
    bin_counts_Data_B_safe(bin_counts_Data_B_safe == 0) = eps;
    
    P_final = bin_counts(2,:) ./ bin_counts_Data_B_safe;
    
    %% Step 2.8: Store Results for Current Iteration
    data_nd_618{iteration_idx, 1} = P_final;       % Store P_final (ratio of counts)
    data_nd_618{iteration_idx, 2} = bin_thresholds;% Store bin thresholds
    
    fprintf('Iteration %d: P_final and bin thresholds stored\n', iteration_idx);
end

%% ===================== Part 3: Save Intermediate Results =====================
save('data_nd_6182_40area_91_40_6_30.mat', 'data_nd_618');


%% ===================== Part 4: Calculate Pdif (Symmetric Bin Differences) =====================
% Reload intermediate results (consistent with original code flow)
load('data_nd_6182_40area2_91_40_6_30.mat');

for iteration_idx = 1:total_iterations
    % Extract P_final for current iteration
    P_final = data_nd_618{iteration_idx, 1};
    half_bin_count = length(P_final) / 2;
    
    % Initialize Pdif matrix (3 rows: P_half+1-i, P_half+i, difference)
    Pdif_matrix = zeros(3, half_bin_count);
    
    for sym_idx = 1:half_bin_count
        % Symmetric bin indices (mirror around half_bin_count)
        left_bin_idx = half_bin_count + 1 - sym_idx;
        right_bin_idx = half_bin_count + sym_idx;
        
        % Extract values and calculate difference
        Pdif_matrix(1, sym_idx) = P_final(left_bin_idx);
        Pdif_matrix(2, sym_idx) = P_final(right_bin_idx);
        Pdif_matrix(3, sym_idx) = Pdif_matrix(1, sym_idx) - Pdif_matrix(2, sym_idx);
    end
    
    % Store Pdif matrix for current iteration
    Pdif{iteration_idx, 1} = Pdif_matrix;
    fprintf('Iteration %d: Pdif calculated (symmetric bins: %d)\n', iteration_idx, half_bin_count);
end

%% ===================== Part 5: Final Save & Completion =====================
% Save final Pdif results (optional: combine with intermediate results)
save('data_nd_6182_40area_91_40_6_30_with_Pdif.mat', 'data_nd_618', 'Pdif');
