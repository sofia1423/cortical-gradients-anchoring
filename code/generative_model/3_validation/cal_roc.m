% =========================================================================
% @function cal_roc
% @brief   Calculate ROC curve metrics (true positive/false positive rates) for neuron weight prediction
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]

% Key Logic:
%   1. Identify zero-degree columns (no connections) in target/ground truth matrices
%   2. Generate ranked list of predicted edges (descending by predicted weight)
%   3. Iteratively threshold predicted edges (1 to 8190) and calculate TP/FP rates
%   4. Filter out zero-degree columns to match model constraints
% =========================================================================

function [true_rate, false_rate, sub_edges] = cal_roc(predicted_log_weights, target_matrix, bin_data, ground_truth_matrix, num_valid_areas)
    %% ===================== Input Validation =====================
    % Check input dimensions to prevent runtime errors
    if size(predicted_log_weights,1) ~= 8190 && size(predicted_log_weights,1) ~= size([bin_data{1};bin_data{2};bin_data{3};bin_data{4}],1)
        error('cal_roc:InvalidInput', 'Predicted log weights must match total edge count (8190)');
    end
    if size(target_matrix) ~= [91,91] || size(ground_truth_matrix) ~= [91,91]
        error('cal_roc:InvalidMatrixSize', 'Target/ground truth matrices must be 91×91');
    end
    
    %% ===================== Initialization =====================
    total_areas = 91;  % Total brain areas in the dataset
    true_rate = zeros(8190, 1);   % True Positive Rate (Sensitivity)
    false_rate = zeros(8190, 1);  % False Positive Rate (1-Specificity)
    sub_edges = zeros(8190, 1);   % Filtered edge count (zero-degree columns removed)
    
    %% ===================== Zero-Degree Column Identification =====================
    % Calculate column degrees for target matrix (m2)
    target_degree = zeros(total_areas, 2);
    target_degree(:,1) = sum(target_matrix')';  % Row degree (sum of rows)
    target_degree(:,2) = sum(target_matrix)';    % Column degree (sum of columns)
    [zero_degree_cols_target, ~] = find(target_degree(:,2) == 0);  % Columns with no connections
    
    % Calculate column degrees for ground truth matrix (m29)
    ground_truth_degree = zeros(total_areas, 2);
    ground_truth_degree(:,1) = sum(ground_truth_matrix')';  % Row degree
    ground_truth_degree(:,2) = sum(ground_truth_matrix)';    % Column degree
    [zero_degree_cols_gt, ~] = find(ground_truth_degree(:,2) == 0);  % Zero-degree columns in ground truth
    
    %% ===================== Extract Target Edges (Non-Zero in m2) =====================
    target_edges = [];
    edge_count = 0;
    
    for source_idx = 1:total_areas
        for target_idx = 1:total_areas
            if target_matrix(source_idx, target_idx) ~= 0
                edge_count = edge_count + 1;
                target_edges(edge_count, 1) = source_idx;    % Source area index
                target_edges(edge_count, 2) = target_idx;    % Target area index
                target_edges(edge_count, 3) = target_matrix(source_idx, target_idx);  % Weight value
            end
        end
    end
    
    %% ===================== Build Full Ranked Edge List =====================
    % Combine edges from all 4 bins (full edge set: 8190 edges)
    full_edges = zeros(8190, 3);
    full_edges(:,1) = [bin_data{1}(:,1); bin_data{2}(:,1); bin_data{3}(:,1); bin_data{4}(:,1)];  % Source indices
    full_edges(:,2) = [bin_data{1}(:,2); bin_data{2}(:,2); bin_data{3}(:,2); bin_data{4}(:,2)];  % Target indices
    full_edges(:,3) = predicted_log_weights;  % Predicted log weights
    
    % Sort edges in DESCENDING order of predicted weight (highest first)
    [~, sorted_indices] = sort(full_edges(:,3), 'descend');
    ranked_edges = full_edges(sorted_indices, :);
    
    %% ===================== Iterative Thresholding (ROC Calculation) =====================
    % Precompute key metrics for ROC calculation
    total_positive = size(target_edges, 1);  % Total actual positive edges (FN + TP)
    total_negative = (total_areas * num_valid_areas) - num_valid_areas - total_positive;  % Total actual negative edges (TN + FP)
    
    % Iterate over all possible edge thresholds (1 to 8190 edges)
    for threshold_idx = 1:8190
        % Extract top N predicted edges (N = threshold_idx)
        top_predicted_edges = ranked_edges(1:threshold_idx, :);
        
        % Initialize adjacency matrix for predicted edges
        predicted_adj_matrix = zeros(total_areas);
        for edge_idx = 1:size(top_predicted_edges, 1)
            src = top_predicted_edges(edge_idx, 1);
            tgt = top_predicted_edges(edge_idx, 2);
            predicted_adj_matrix(src, tgt) = top_predicted_edges(edge_idx, 3);
        end
        
        %% Filter predicted matrix (remove zero-degree columns from target matrix)
        filtered_pred_matrix = predicted_adj_matrix;
        for col_idx = 1:size(zero_degree_cols_target, 1)
            filtered_pred_matrix(:, zero_degree_cols_target(col_idx)) = 0;
        end
        
        %% Filter predicted matrix (remove zero-degree columns from ground truth)
        filtered_pred_matrix_gt = predicted_adj_matrix;
        for col_idx = 1:size(zero_degree_cols_gt, 1)
            filtered_pred_matrix_gt(:, zero_degree_cols_gt(col_idx)) = 0;
        end
        
        %% Calculate TP/FP/TN/FN metrics
        true_positives = 0;    % TP: Predicted edge exists in target matrix
        total_predicted = 0;   % Total predicted edges (TP + FP)
        total_predicted_gt = 0;% Total predicted edges (ground truth filtered)
        
        for source_idx = 1:total_areas
            for target_idx = 1:total_areas
                % Count predicted edges (ground truth filtered)
                if filtered_pred_matrix_gt(source_idx, target_idx) ~= 0
                    total_predicted_gt = total_predicted_gt + 1;
                end
                
                % Count predicted edges (target matrix filtered)
                if filtered_pred_matrix(source_idx, target_idx) ~= 0
                    total_predicted = total_predicted + 1;
                    % Check if predicted edge is a true positive
                    if target_matrix(source_idx, target_idx) ~= 0
                        true_positives = true_positives + 1;
                    end
                end
            end
        end
        
        %% Calculate ROC metrics
        % True Positive Rate (Sensitivity) = TP / (TP + FN)
        if total_positive > 0
            true_rate(threshold_idx, 1) = true_positives / total_positive;
        else
            true_rate(threshold_idx, 1) = 0;
        end
        
        % False Positive Rate (1 - Specificity) = FP / (FP + TN)
        false_positives = total_predicted - true_positives;
        if total_negative > 0
            false_rate(threshold_idx, 1) = (total_negative - (total_predicted - true_positives)) / total_negative;
        else
            false_rate(threshold_idx, 1) = 0;
        end
        
        % Store filtered edge count (ground truth zero-degree columns removed)
        sub_edges(threshold_idx, 1) = total_predicted_gt;
    end
    
    %% ===================== Final Sanity Checks =====================
    % Ensure ROC metrics are within valid range [0,1]
    true_rate(true_rate < 0) = 0;
    true_rate(true_rate > 1) = 1;
    false_rate(false_rate < 0) = 0;
    false_rate(false_rate > 1) = 1;
    
end

