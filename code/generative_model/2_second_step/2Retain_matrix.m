% =========================================================================
% @file    neuron_correlation_analysis_40area.m
% @brief   Neuron Connectivity Correlation Analysis (40-Area) with Random Sampling
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key workflow:
%          1. Load core datasets (gradient, distance, cluster, weight matrices)
%          2. Revise P_final with scaling factors (1.05-2.0)
%          3. Generate 24-category p_sln2 matrix from gradient ratios
%          4. Random sampling (1000 iterations) with weight matrix generation
%          5. Calculate correlations (cor2-cor7) and filter results (cor2>0.43)
%          6. Save filtered results (rangepre + full correlation data)
% =========================================================================

%% ===================== Part 1: Environment Setup & Data Loading =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;
rng('shuffle');  % Set random seed (consistent random sampling)

% Define core parameters (centralized for easy modification)
mm27 = 1000;                % Number of random sampling iterations
bin_count = 24;             % Number of bins for category_24 function
kkk_start = 2;              % Start index for kkk loop
kkk_end = 20;               % End index for kkk loop
cor2_threshold = 0.43;      % Core filtering threshold (cor2 > 0.43)

% Load core datasets (critical dependencies)

load('data_nd_6182_40area_91_40_6_30.mat');  % P_final + thr_Bin data
load('distance_matrix.mat');                 % Pairwise distance matrix (dist_m)
load('axis3_cluster.mat');                   % Cluster data (axis3_cluster)
load('lp19_gd_fitnew.mat');                  % Gradient fit results (gdfit10)
load('data4_F_40area_40_40.mat');            % 40-area weight feature data
load('MatF_SLN2021.mat');  % 91x91 weight matrix (MatF)
load('efferent_afferent_newp40_replace40.mat');  % Efferent/afferent data (range_newp)


%% ===================== Part 2: Preprocess Core Data =====================
% Extract P_final and threshold data for iteration 10
P_final1 = data_nd_618{10,1};
thr_Bin_final = data_nd_618{10,2};

% Revise P_final with scaling factors (1.05 to 2.0, step 0.05: 20 values)
fprintf('Generating revised P_final values (scaling 1.05-2.0)...\n');
P_revise66 = cell(20,1);  % Pre-allocate cell array
for scale_idx = 1:20
    scaling_factor = 1 + 0.05 * scale_idx;  % 1.05, 1.10, ..., 2.00
    P_revise66{scale_idx} = P_final1 * scaling_factor;
end

% Compute 91x91 gradient ratio matrix (log(gdfit10(k,1)/gdfit10(i,1)))
fprintf('Computing 91x91 gradient ratio matrix...\n');
gradient_ratio_matrix = zeros(91,91);
for source_area = 1:91
    for target_area = 1:91
        gradient_ratio_matrix(source_area, target_area) = log(gdfit10(target_area, 1) / gdfit10(source_area, 1));
    end
end
m22 = gradient_ratio_matrix;  % Retain original variable name for compatibility

% Preprocess weight feature data (remove columns 3-9)
weight_feature_data = data4_F;
weight_feature_data(:,3:9) = [];

% Sort weight feature data by column 3 (descending) using custom matrix_sort function
fprintf('Sorting weight feature data...\n');
weight_r = matrix_sort(weight_feature_data, 3, 1, 'descend');

% Rename 91x91 weight matrix for clarity
weight91_n = MatF;
range_newp = efferent_afferent;


%% ===================== Part 3: Random Sampling & Correlation Calculation =====================
% Pre-allocate storage for results (critical for performance)
data_f = cell(kkk_end, mm27);
data_f2 = cell(kkk_end, mm27);
matrix_k = cell(kkk_end, mm27);
matrix_k_nor = cell(mm27,1);
range_pre = cell(kkk_end, mm27);

% Correlation storage (pre-allocate matrices for speed)
cor2 = zeros(kkk_end, mm27); p2 = zeros(kkk_end, mm27);
cor3 = zeros(kkk_end, mm27); p3 = zeros(kkk_end, mm27);
cor4 = zeros(kkk_end, mm27); p4 = zeros(kkk_end, mm27);
cor5 = zeros(kkk_end, mm27); p5 = zeros(kkk_end, mm27);
cor6 = zeros(kkk_end, mm27); p6 = zeros(kkk_end, mm27);
cor7 = zeros(kkk_end, mm27); p7 = zeros(kkk_end, mm27);

% Main loop over revised P_final values (kkk=2 to 20)
fprintf('Starting random sampling & correlation analysis (total iterations: %d x %d = %d)\n', ...
    kkk_end-kkk_start+1, mm27, (kkk_end-kkk_start+1)*mm27);

progress_counter = 0;
total_iterations = (kkk_end - kkk_start + 1) * mm27;

for kkk = kkk_start:kkk_end
    % Get revised P_final for current kkk
    P_revise = P_revise66{kkk};
    
    % Generate 24-category p_sln2 matrix
    p_sln2 = category_24(m22, thr_Bin_final, P_revise(1,:));
    
    % Assign p_sln2 values to weight_r column 4
    for row_idx = 1:8190
        weight_r(row_idx,4) = p_sln2(weight_r(row_idx,1), weight_r(row_idx,2));
    end
    
    % Parallel random sampling (parfor for speed)
    parfor sampling_idx = 1:mm27
        [data_f{kkk,sampling_idx}, data_f2{kkk,sampling_idx}] = rand_retain_train_40_40(weight_r, weight91_n);
    end
    
    % Process sampling results (calculate correlations)
    for sampling_idx = 1:mm27
        progress_counter = progress_counter + 1;
        
        % Only process valid samples (data_f2 > 0)
        if data_f2{kkk,sampling_idx} > 0
            % Generate link matrix (91x91)
            matrix_k{kkk,sampling_idx} = link_matrix(data_f{kkk,sampling_idx}(1:data_f2{kkk,sampling_idx},1:3), 91);
            
            % Normalize link matrix
            matrix_k_nor{sampling_idx} = nor_matrix(matrix_k{kkk,sampling_idx});
            
            % Calculate range for 40 areas
            [range_pre{kkk,sampling_idx}(:,1), range_pre{kkk,sampling_idx}(:,2)] = cal_range_40area(matrix_k_nor{sampling_idx}, dist_m);
            
            % Calculate correlations (cor2-cor7) with p-values
            [cor2(kkk,sampling_idx), p2(kkk,sampling_idx)] = corr_lp(range_pre{kkk,sampling_idx}(:,2), range_newp(:,2));
            [cor3(kkk,sampling_idx), p3(kkk,sampling_idx)] = corr_lp(range_pre{kkk,sampling_idx}(:,2), axis3_cluster(:,1));
            [cor4(kkk,sampling_idx), p4(kkk,sampling_idx)] = corr_lp(range_pre{kkk,sampling_idx}(:,2), axis3_cluster(:,2));
            [cor5(kkk,sampling_idx), p5(kkk,sampling_idx)] = corr_lp(range_pre{kkk,sampling_idx}(:,2), axis3_cluster(:,3));
            [cor6(kkk,sampling_idx), p6(kkk,sampling_idx)] = corr_lp(range_pre{kkk,sampling_idx}(:,2), axis3_cluster(:,4));
            [cor7(kkk,sampling_idx), p7(kkk,sampling_idx)] = corr_lp(range_pre{kkk,sampling_idx}(:,2), axis3_cluster(:,5));
        end
        
        % Progress update every 100 iterations
        if mod(progress_counter, 100) == 0
            fprintf('Progress: %d/%d iterations completed (%.1f%%)\n', ...
                progress_counter, total_iterations, (progress_counter/total_iterations)*100);
        end
    end
end


%% ===================== Part 4: Filter & Save Results (cor2 > 0.43) =====================
% Initialize storage for filtered results
result_save = {};  % Cell array: store full results (struct)
rangepre_all = []; % Matrix: collect all valid rangepre values
valid_result_count = 0;  % Counter for valid results

fprintf('Filtering results (cor2 > %.2f)...\n', cor2_threshold);

% Double loop to filter valid results
for kkk = kkk_start:kkk_end
    for sampling_idx = 1:mm27
        % Only process valid samples (data_f2 > 0)
        if data_f2{kkk,sampling_idx} > 0
            % Core filtering condition: cor2 > threshold
            if cor2(kkk,sampling_idx) > cor2_threshold
                valid_result_count = valid_result_count + 1;
                
                % Extract all key data for current valid result
                current_rangepre = range_pre{kkk,sampling_idx}(:,2);
                current_cor = [cor2(kkk,sampling_idx), cor3(kkk,sampling_idx), cor4(kkk,sampling_idx), ...
                               cor5(kkk,sampling_idx), cor6(kkk,sampling_idx), cor7(kkk,sampling_idx)];
                current_p = [p2(kkk,sampling_idx), p3(kkk,sampling_idx), p4(kkk,sampling_idx), ...
                             p5(kkk,sampling_idx), p6(kkk,sampling_idx), p7(kkk,sampling_idx)];
                current_data_f = data_f{kkk,sampling_idx};
                current_data_f2 = data_f2{kkk,sampling_idx};
                current_index = [kkk, sampling_idx];  % Backtracking index
                
                % Store in struct (clear naming for downstream analysis)
                result_save{valid_result_count} = struct(...
                    'cor_values', current_cor, ...
                    'p_values', current_p, ...
                    'rangepre', current_rangepre, ...
                    'data_f', current_data_f, ...
                    'data_f2', current_data_f2, ...
                    'index_kkk_sampling', current_index);
                
                % Append rangepre to collective matrix (row-wise)
                rangepre_all = [rangepre_all; current_rangepre'];
            end
        end
    end
end

% Save results (only if valid results exist)
if valid_result_count > 0
    % Save filtered results
    save('rangepre_40area_40_40weight_91_40P.mat', 'rangepre_all');
    save('matrix_save_40area_40_40weight_91_40P.mat', 'result_save', 'valid_result_count');
    
    % Print success message
    fprintf('✅ Analysis completed! %d valid results found (cor2 > %.2f)\n', valid_result_count, cor2_threshold);
    fprintf('Results saved to:\n');
    fprintf('  - rangepre_40area_40_40weight_91_40P.mat (rangepre_all)\n');
    fprintf('  - matrix_save_40area_40_40weight_91_40P.mat (result_save + count)\n');
else
    % Print no results message
    fprintf('❌ Analysis completed! No valid results found (cor2 > %.2f)\n', cor2_threshold);
end
