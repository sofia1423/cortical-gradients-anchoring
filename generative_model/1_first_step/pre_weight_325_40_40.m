% =========================================================================
% Neuron Connectivity Analysis: 2x2 Binning & Linear Fitting
%% 1. Load Data & Preprocess Weight Matrix
% Load core datasets (neuron connectivity/gradient/distance data)
load('MatF_SLN2021.mat');          % Neuron weight matrix (MatF)
load('lp19_gd_fitnew.mat');        % Gradient fit results (gdfit10)
load('distance_matrix.mat');       % Pairwise distance matrix (dist_m)

% Remove zero columns from weight matrix (filter valid 91-area indices)
zero_cols = find(all(MatF == 0, 1));       % Column indices with all zero values
area_indices_91 = 1:91;                    % Initial full 91-area index range
area_indices_91(zero_cols) = [];           % Remove zero columns (valid areas only)
weight_matrix_40x40 = trans_matrix_40area(MatF, 40, 40);  % Transform to 40x40 weight matrix

%% 2. Build Feature Matrix (non-zero weight entries with gradient/distance)
% Initialize feature matrix to store core neuron connectivity features
clear neuron_feature_matrix;
feature_row_count = 0;

% Iterate through non-zero entries in weight matrix
for weight_row = 1:size(weight_matrix_40x40, 1)
    for weight_col = 1:size(weight_matrix_40x40, 2)
        if weight_matrix_40x40(weight_row, weight_col) ~= 0
            feature_row_count = feature_row_count + 1;
            
            % 1. Core indices (original weight matrix + mapped 40-area indices)
            neuron_feature_matrix(feature_row_count, 1)  = weight_row;                  % Row index in weight matrix
            neuron_feature_matrix(feature_row_count, 2)  = weight_col;                  % Column index in weight matrix
            neuron_feature_matrix(feature_row_count, 10) = area_indices_91(weight_row); % Mapped row index (40-area)
            neuron_feature_matrix(feature_row_count, 11) = area_indices_91(weight_col); % Mapped column index (40-area)
            
            % 2. Weight-related features
            neuron_feature_matrix(feature_row_count, 3)  = log(weight_matrix_40x40(weight_row, weight_col)); % Log-transformed weight
            neuron_feature_matrix(feature_row_count, 4)  = weight_matrix_40x40(weight_row, weight_col);       % Raw weight value
            
            % 3. Gradient & distance features
            grad_ratio = log(gdfit10(area_indices_91(weight_row), 1) / gdfit10(area_indices_91(weight_col), 1));
            neuron_feature_matrix(feature_row_count, 5)  = abs(grad_ratio);             % Absolute log gradient ratio
            neuron_feature_matrix(feature_row_count, 6)  = dist_m(area_indices_91(weight_row), area_indices_91(weight_col)); % Pairwise distance
        end
    end
end

%% 3. 2x2 Binning & Linear Fitting Pipeline
% Initialize output matrix to store binning/fitting results
clear fitting_results; 
valid_bin_count = 0;

% Outer loop: Distance threshold sweep (0.1 to 70, step 0.1)
for dist_thresh_idx = 1:700
    distance_threshold = 0.1 * dist_thresh_idx;  % Distance split value (column 6 of feature matrix)
    
    % Split feature matrix by distance threshold (logical indexing for speed)
    is_low_distance = neuron_feature_matrix(:,6) <= distance_threshold;
    is_high_distance = neuron_feature_matrix(:,6) > distance_threshold;
    feature_low_distance = neuron_feature_matrix(is_low_distance, :);
    feature_high_distance = neuron_feature_matrix(is_high_distance, :);

    % Inner loop: Gradient threshold sweep (0.01 to 1, step 0.01)
    for grad_thresh_idx = 1:100
        gradient_threshold = 0.01 + 0.01 * (grad_thresh_idx - 1);  % Gradient split value (column 5)
        
        % Split low/high distance data by gradient threshold
        is_low_grad_low_dist = feature_low_distance(:,5) <= gradient_threshold;
        is_high_grad_low_dist = feature_low_distance(:,5) > gradient_threshold;
        is_low_grad_high_dist = feature_high_distance(:,5) <= gradient_threshold;
        is_high_grad_high_dist = feature_high_distance(:,5) > gradient_threshold;
        
        % 2x2 binning (4 bins: distance × gradient combinations)
        bin_low_dist_low_grad = feature_low_distance(is_low_grad_low_dist, :);   % Bin 1: Low dist + Low grad
        bin_low_dist_high_grad = feature_low_distance(is_high_grad_low_dist, :);% Bin 2: Low dist + High grad
        bin_high_dist_low_grad = feature_high_distance(is_low_grad_high_dist, :);% Bin3: High dist + Low grad
        bin_high_dist_high_grad = feature_high_distance(is_high_grad_high_dist, :);% Bin4: High dist + High grad

        % Only process if all bins have >1 sample (statistical validity)
        bin1_size = size(bin_low_dist_low_grad, 1);
        bin2_size = size(bin_low_dist_high_grad, 1);
        bin3_size = size(bin_high_dist_low_grad, 1);
        bin4_size = size(bin_high_dist_high_grad, 1);
        
        if (bin1_size > 1) && (bin2_size > 1) && (bin3_size > 1) && (bin4_size > 1)
            valid_bin_count = valid_bin_count + 1;
            
            % Store bin metadata (thresholds + sample counts)
            fitting_results(valid_bin_count, 11) = distance_threshold;  % Distance split threshold
            fitting_results(valid_bin_count, 12) = gradient_threshold;  % Gradient split threshold
            fitting_results(valid_bin_count, 1:4) = [bin1_size, bin2_size, bin3_size, bin4_size]; % Sample count per bin

            % Initialize arrays for global correlation calculation
            fitted_log_weights = [];  % Fitted log-weight values from regression
            original_log_weights = [];% Original log-weight values from feature matrix

            % Process each bin (normalization + linear regression)
            bin_list = {bin_low_dist_low_grad, bin_low_dist_high_grad, bin_high_dist_low_grad, bin_high_dist_high_grad};
            for bin_index = 1:4
                bin_data = bin_list{bin_index};
                
                % Min-Max normalization (scale features to [0,1] for regression)
                bin_data(:,8) = (bin_data(:,5) - min(bin_data(:,5))) / (max(bin_data(:,5)) - min(bin_data(:,5))); % Normalized gradient
                bin_data(:,9) = (bin_data(:,6) - min(bin_data(:,6))) / (max(bin_data(:,6)) - min(bin_data(:,6))); % Normalized distance

                % Linear regression: log(weight) ~ normalized(gradient + distance)
                regression_model = fitlm(bin_data(:,8:9), bin_data(:,3));

                % Store regression performance metrics (sqrt of R² for interpretability)
                fitting_results(valid_bin_count, bin_index + 4) = sqrt(regression_model.Rsquared.Ordinary);
                
                % Extract & store regression coefficients (gradient + distance)
                grad_coefficient = table2array(regression_model.Coefficients(2,1));  % Gradient feature coefficient
                dist_coefficient = table2array(regression_model.Coefficients(3,1));  % Distance feature coefficient
                fitting_results(valid_bin_count, bin_index + (bin_index-1)*2 + 12) = grad_coefficient;
                fitting_results(valid_bin_count, bin_index + (bin_index-1)*2 + 13) = dist_coefficient;

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
                fitting_results(valid_bin_count, 9:10) = NaN;  % Handle empty data case (avoid errors)
            end
        end
    end
end

%% 4. Filter Results (Negative Coefficients Only)
% Select rows where ALL 8 coefficients (4 bins × 2 features) are negative
% Column indices for coefficients: 13(grad1),14(dist1),16(grad2),17(dist2),19(grad3),20(dist3),22(grad4),23(dist4)
is_negative_coeff = ...
    (fitting_results(:,13)<0) & (fitting_results(:,14)<0) & ...
    (fitting_results(:,16)<0) & (fitting_results(:,17)<0) & ...
    (fitting_results(:,19)<0) & (fitting_results(:,20)<0) & ...
    (fitting_results(:,22)<0) & (fitting_results(:,23)<0);

negative_coeff_indices = find(is_negative_coeff);
final_results_negative_coeff = fitting_results(negative_coeff_indices, :);  % Final filtered results