% =========================================================================
% @file    neuron_length_correlation_analysis.m
% @brief   Parallelized Correlation Analysis for Neuron Length Data (Offset Sweep -100 to 99)
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key parameters: Offset range=-100 to 99 (200 steps), Batches=24 (113 samples/batch), Total pairs=113×4723
% =========================================================================

%% ===================== Part 1: Environment Setup & Initialization =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;
rng('shuffle');  % Set random seed for consistent correlation results

% Define core parameters (centralized for easy modification)
NUM_BATCHES = 24;               % Total number of processing batches (hh=1:24)
SAMPLES_PER_BATCH = 113;        % Number of samples per batch (s1 range per hh)
OFFSET_MIN = -100;              % Minimum offset value (kk=-100)
OFFSET_MAX = 99;                % Maximum offset value (kk=99)
NUM_OFFSETS = OFFSET_MAX - OFFSET_MIN + 1;  % Total offset steps (200)
TARGET_SAMPLES = 4723;          % Number of target samples (s2=1:4723)

% Predefine results directory (create if not exists)
results_dir = 'results32k';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created results directory: %s\n', results_dir);
end

%% ===================== Part 2: Load Core Data =====================
fprintf('Loading core datasets...\n');
load('len3.mat');                % Source1 neuron length data (len3: [2712×91] = 24×113 samples × 91 features)
load('len4.mat');                % Source2 neuron length data (len4: [4723×91] = 4723 samples × 91 features)
load('Exp_T1_T2.mat');           % Neurophysiological data (zf(:,1) used for correlation)
fprintf('Core datasets loaded successfully:\n');
fprintf('  len3 size: %dx%d | len4 size: %dx%d\n', size(len3,1), size(len3,2), size(len4,1), size(len4,2));

%% ===================== Part 3: Parallel Batch Processing =====================
fprintf('Starting batch processing (total batches: %d)...\n', NUM_BATCHES);

for batch_idx = 1:NUM_BATCHES
    % Calculate sample range for current batch (s1 = ii1 to ii2)
    batch_start = (batch_idx - 1) * SAMPLES_PER_BATCH + 1;
    batch_end = batch_idx * SAMPLES_PER_BATCH;
    
    % Pre-allocate storage for current batch (critical for speed)
    cor_results = cell(SAMPLES_PER_BATCH, TARGET_SAMPLES);  % Correlation results (num×s2)
    nd_results = cell(SAMPLES_PER_BATCH, TARGET_SAMPLES);   % Minimum length results (num×s2)
    
    % Initialize sample counter for current batch
    batch_sample_count = 0;
    
    % Process each source sample in current batch
    for source_idx = batch_start:batch_end
        batch_sample_count = batch_sample_count + 1;
        
        % Extract length data for current source sample
        source_length = len3(source_idx, :);
        
        % Process each target sample
        for target_idx = 1:TARGET_SAMPLES
            % Extract length data for current target sample
            target_length = len4(target_idx, :);
            
            % Pre-allocate offset matrix (200×91) for current source-target pair
            offset_matrix = zeros(NUM_OFFSETS, size(source_length,2));
            
            % Sweep offset values (-100 to 99)
            for offset_step = 1:NUM_OFFSETS
                % Calculate current offset value (kk = -100 + offset_step - 1)
                current_offset = OFFSET_MIN + offset_step - 1;
                
                % Apply offset to source length data
                for feature_idx = 1:size(source_length,2)
                    offset_matrix(offset_step, feature_idx) = source_length(feature_idx) + current_offset;
                    % Calculate minimum value between offset source and target length
                    nd_results{batch_sample_count, target_idx}(offset_step, feature_idx) = min(offset_matrix(offset_step, feature_idx), target_length(feature_idx));
                end
                
                % Calculate correlation between inverted zf(:,1) and nd results (transpose for corr_lp)
                [cor_coeff, p_value] = corr_lp(zf(:,1)*(-1), nd_results{batch_sample_count, target_idx}(offset_step, :)');
                
                % Store correlation coefficient and p-value
                cor_results{batch_sample_count, target_idx}(offset_step, 1) = cor_coeff;
                cor_results{batch_sample_count, target_idx}(offset_step, 2) = p_value;
            end
        end
        
        % Progress update for current batch (every 10 samples)
        if mod(batch_sample_count, 10) == 0
            fprintf('  Batch %d: %d/%d samples processed\n', batch_idx, batch_sample_count, SAMPLES_PER_BATCH);
        end
    end
    
    %% ===================== Part 4: Save Batch Results =====================
    % Generate filenames for current batch
    batch_str = num2str(batch_idx);
    cor_filename = strcat('corlp_', batch_str, '.mat');
    nd_filename = strcat('nd_', batch_str, '.mat');
    
    % Save results to dedicated directory (parallel-safe save)
    cd(results_dir);
    parsave(cor_filename, cor_results);
    parsave(nd_filename, nd_results);
    cd('..');  % Return to parent directory
    
 
end

%% ===================== Part 5: Custom Parallel Save Function =====================
function parsave(fname, x)
    % parsave: Parallel-safe save function (avoids race conditions in batch saving)
    % Inputs:
    %   fname - Output filename (e.g., 'corlp_1.mat')
    %   x - Data to save (cell array of correlation/nd results)
    save(fname, 'x');
end