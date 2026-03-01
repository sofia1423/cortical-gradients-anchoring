% =========================================================================
% Neuron Edge Weight Prediction (8190 Edges)
% Purpose: Predict full 91x91 neuron connectivity weights using 4-bin linear fitting rules
%          derived from non-zero weight entries, then validate predictions
% =========================================================================

%% 1. Load & Preprocess Core Data
% Load neuron weight matrix, gradient fit, and distance matrix
load('MatF_SLN2021.mat');
load('lp19_gd_fitnew.mat');          % Gradient fit data (gdfit10)
load('distance_matrix.mat');         % Pairwise distance matrix (dist_m)

% Remove zero columns from weight matrix (filter valid 91-area indices)
zero_columns = find(all(MatF == 0, 1));
area_indices_91 = 1:91;
area_indices_91(zero_columns) = [];
weight_matrix_40x40 = trans_matrix_40area(MatF, 40, 40);  % Transform to 40x40 weight matrix

%% 2. Build Feature Matrix for Non-Zero Weights (Training Data)
clear neuron_feature_matrix;
feature_row_count = 0;

% Extract non-zero weight entries with gradient/distance features
for weight_row = 1:size(weight_matrix_40x40, 1)
    for weight_col = 1:size(weight_matrix_40x40, 2)
        if weight_matrix_40x40(weight_row, weight_col) ~= 0
            feature_row_count = feature_row_count + 1;
            
            % Core indices (original + mapped 91-area indices)
            neuron_feature_matrix(feature_row_count, 1)  = weight_row;
            neuron_feature_matrix(feature_row_count, 2)  = weight_col;
            neuron_feature_matrix(feature_row_count, 10) = area_indices_91(weight_row);
            neuron_feature_matrix(feature_row_count, 11) = area_indices_91(weight_col);
            
            % Weight features (log-transformed + raw)
            neuron_feature_matrix(feature_row_count, 3)  = log(weight_matrix_40x40(weight_row, weight_col));
            neuron_feature_matrix(feature_row_count, 4)  = weight_matrix_40x40(weight_row, weight_col);
            
            % Gradient & distance features
            grad_ratio = log(gdfit10(area_indices_91(weight_row), 1) / gdfit10(area_indices_91(weight_col), 1));
            neuron_feature_matrix(feature_row_count, 5)  = abs(grad_ratio);  % Absolute log gradient ratio
            neuron_feature_matrix(feature_row_count, 6)  = dist_m(area_indices_91(weight_row), area_indices_91(weight_col));
        end
    end
end

%% 3. Define 4 Bins (Distance + Gradient Thresholds for 40x40 Version)
% Thresholds (2026.2.7 update: 40area x 40area)
distance_threshold = 33.5;
gradient_threshold = 0.19;

% Bin indices (split by distance + gradient thresholds)
bin_indices{1} = find((neuron_feature_matrix(:,6) <= distance_threshold) & (neuron_feature_matrix(:,5) <= gradient_threshold));
bin_indices{2} = find((neuron_feature_matrix(:,6) > distance_threshold) & (neuron_feature_matrix(:,5) <= gradient_threshold));
bin_indices{3} = find((neuron_feature_matrix(:,6) <= distance_threshold) & (neuron_feature_matrix(:,5) > gradient_threshold));
bin_indices{4} = find((neuron_feature_matrix(:,6) > distance_threshold) & (neuron_feature_matrix(:,5) > gradient_threshold));

%% 4. Train Linear Models for Each Bin (Normalization + Fitting)
clear bin_data_list;
clear regression_models;
clear normalization_params;  % Store min/max for normalization (per bin)
fitted_log_weights = [];
original_log_weights = [];

for bin_idx = 1:4
    % Extract data for current bin
    bin_data_list{bin_idx} = neuron_feature_matrix(bin_indices{bin_idx}, :);
    
    % Min-Max normalization for gradient (column 5)
    grad_min = min(bin_data_list{bin_idx}(:,5));
    grad_range = max(bin_data_list{bin_idx}(:,5)) - grad_min;
    bin_data_list{bin_idx}(:,8) = (bin_data_list{bin_idx}(:,5) - grad_min) / grad_range;
    normalization_params(bin_idx, 1:2) = [grad_min, grad_range];  % Store gradient norm params
    
    % Min-Max normalization for distance (column 6)
    dist_min = min(bin_data_list{bin_idx}(:,6));
    dist_range = max(bin_data_list{bin_idx}(:,6)) - dist_min;
    bin_data_list{bin_idx}(:,9) = (bin_data_list{bin_idx}(:,6) - dist_min) / dist_range;
    normalization_params(bin_idx, 3:4) = [dist_min, dist_range];  % Store distance norm params
    
    % Train linear regression model (log weight ~ normalized gradient + distance)
    regression_models{bin_idx} = fitlm(bin_data_list{bin_idx}(:,8:9), bin_data_list{bin_idx}(:,3));
    
    % Accumulate fitted/original values for correlation calculation
    fitted_log_weights = [fitted_log_weights; regression_models{bin_idx}.Fitted];
    original_log_weights = [original_log_weights; bin_data_list{bin_idx}(:,3)];
end

% Calculate correlation between fitted and original log weights (training data)
[train_corr, train_p_value] = corr_lp(original_log_weights, fitted_log_weights);
fprintf('Training correlation (non-zero weights): %.4f\n', train_corr);  % Expected ~0.5329

%% 5. Predict Weights for All 8190 Edges (Full 91x91 Matrix)
clear full_edge_matrix;
edge_count = 0;

% Build matrix for all non-self edges (91x91, exclude i=k)
for area_i = 1:91
    for area_k = 1:91
        if area_i ~= area_k
            edge_count = edge_count + 1;
            full_edge_matrix(edge_count, 1) = area_i;          % Source area
            full_edge_matrix(edge_count, 2) = area_k;          % Target area
            full_edge_matrix(edge_count, 3) = abs(log(gdfit10(area_i,1)/gdfit10(area_k,1)));  % Gradient ratio
            full_edge_matrix(edge_count, 4) = dist_m(area_i, area_k);  % Pairwise distance
        end
    end
end

% Assign full edges to 4 bins (same thresholds as training)
full_bin_indices{1} = find((full_edge_matrix(:,4) <= distance_threshold) & (full_edge_matrix(:,3) <= gradient_threshold));
full_bin_indices{2} = find((full_edge_matrix(:,4) > distance_threshold) & (full_edge_matrix(:,3) <= gradient_threshold));
full_bin_indices{3} = find((full_edge_matrix(:,4) <= distance_threshold) & (full_edge_matrix(:,3) > gradient_threshold));
full_bin_indices{4} = find((full_edge_matrix(:,4) > distance_threshold) & (full_edge_matrix(:,3) > gradient_threshold));

% Load pre-saved normalization parameters (or use in-memory params)
% load('normal40area_40_40.mat');  % Uncomment if using saved params

% Predict log weights for each bin using trained models
clear full_bin_data;
for bin_idx = 1:4
    full_bin_data{bin_idx} = full_edge_matrix(full_bin_indices{bin_idx}, :);
    
    % Apply training-set normalization (critical for consistent prediction)
    full_bin_data{bin_idx}(:,8) = (full_bin_data{bin_idx}(:,3) - normalization_params(bin_idx,1)) / normalization_params(bin_idx,2);
    full_bin_data{bin_idx}(:,9) = (full_bin_data{bin_idx}(:,4) - normalization_params(bin_idx,3)) / normalization_params(bin_idx,4);
    
    % Predict log weights using regression coefficients
    intercept = table2array(regression_models{bin_idx}.Coefficients(1,1));
    grad_coeff = table2array(regression_models{bin_idx}.Coefficients(2,1));
    dist_coeff = table2array(regression_models{bin_idx}.Coefficients(3,1));
    full_bin_data{bin_idx}(:,10) = intercept + full_bin_data{bin_idx}(:,8)*grad_coeff + full_bin_data{bin_idx}(:,9)*dist_coeff;
end

% Combine all predicted edges into one matrix
predicted_full_edges = [full_bin_data{1}; full_bin_data{2}; full_bin_data{3}; full_bin_data{4}];

%% 6. Validate Predictions Against Non-Zero Training Weights
clear prediction_indices;
validation_edge_count = size(neuron_feature_matrix, 1);  % ~999 edges

% Match predicted edges to training edges (by source/target area)
for idx = 1:validation_edge_count
    source_area = neuron_feature_matrix(idx, 10);
    target_area = neuron_feature_matrix(idx, 11);
    prediction_indices(idx,1) = find((predicted_full_edges(:,1)==source_area) & (predicted_full_edges(:,2)==target_area));
end

% Extract predicted log weights for validation
predicted_log_weights = predicted_full_edges(prediction_indices, 10);
actual_log_weights = neuron_feature_matrix(:,3);

% Calculate validation correlation
[validation_corr, validation_p_value] = corr_lp(actual_log_weights, predicted_log_weights);
fprintf('Validation correlation (training edges): %.4f\n', validation_corr);  % Expected ~0.5329

% Z-score normalization for comparison
pred_actual_matrix = [predicted_log_weights, actual_log_weights];
pred_actual_matrix(:,1) = zscore(pred_actual_matrix(:,1));
pred_actual_matrix(:,2) = zscore(pred_actual_matrix(:,2));

%% 7. Build Full 91x91 Predicted Weight Matrix
clear predicted_weight_matrix;

% Assign predicted log weights to 91x91 matrix
for idx = 1:size(predicted_full_edges, 1)
    source_area = predicted_full_edges(idx, 1);
    target_area = predicted_full_edges(idx, 2);
    predicted_weight_matrix(source_area, target_area) = predicted_full_edges(idx, 10);
end

% Save full predicted matrix
save('matrix_pre_8190_40area.mat', 'predicted_weight_matrix');
