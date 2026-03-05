% =========================================================================
% @file    brain_region_shortest_path_analysis.m
% @brief   Parallelized Shortest Path Calculation for 91 Brain Regions (10020 Source Vertices)
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key parameters: Total vertices=10020, Batches=167 (60 vertices/batch), Target regions=91, Results dir=result10020
%          Reference: myelin (83,39), neuron density (64,11)
% =========================================================================

%% ===================== Part 1: Environment Setup & Initialization =====================
% Clear workspace/command window for reproducibility
clearvars; clc; close all;

% Define core parameters (centralized for easy modification)
NUM_BATCHES = 167;               % Total processing batches (i=1:167)
VERTICES_PER_BATCH = 60;         % Vertices processed per batch (60/batch)
TOTAL_VERTICES = NUM_BATCHES * VERTICES_PER_BATCH;  % Total source vertices (10020)
TARGET_REGIONS = 91;             % Number of target brain regions (k=1:91)
RESULTS_DIR = 'result10020';     % Directory for saving batch results

% Create results directory if it does not exist (prevent runtime errors)
if ~exist(RESULTS_DIR, 'dir')
    mkdir(RESULTS_DIR);
    fprintf('Created results directory: %s\n', RESULTS_DIR);
else
    fprintf('Using existing results directory: %s\n', RESULTS_DIR);
end

%% ===================== Part 2: Load Core Data =====================
fprintf('Loading core datasets (graph + vertex indices + region vertices)...\n');
load('G_vertex.mat');                  % Graph object (G) for shortest path calculation
load('indf_10020source.mat');          % Source vertex indices (indf: 10020x1)
load('vertex_of_91regions_140k.mat');  % Vertex IDs for 91 brain regions (min_id: 91x3)
fprintf('Core datasets loaded successfully:\n');
fprintf('  Total source vertices: %d | Target regions: %d\n', TOTAL_VERTICES, TARGET_REGIONS);

%% ===================== Part 3: Parallel Batch Processing (Shortest Path Calculation) =====================
fprintf('Starting shortest path calculation (total batches: %d)...\n', NUM_BATCHES);

% Parallel batch loop (parfor for multi-core acceleration)
parfor batch_idx = 1:NUM_BATCHES
    % Calculate vertex range for current batch (avoid magic numbers)
    batch_start_idx = 1 + VERTICES_PER_BATCH * (batch_idx - 1);
    batch_end_idx = VERTICES_PER_BATCH * batch_idx;
    
    % Pre-allocate shortest path matrix (critical for speed: 60x91 zeros)
    shortest_path_matrix = zeros(VERTICES_PER_BATCH, TARGET_REGIONS);
    
    % Initialize vertex counter for current batch
    batch_vertex_count = 0;
    
    % Process each vertex in current batch
    for vertex_global_idx = batch_start_idx:batch_end_idx
        batch_vertex_count = batch_vertex_count + 1;
        
        % Get source vertex ID from index list
        source_vertex_id = indf(vertex_global_idx, 1);
        
        % Calculate shortest path to each target brain region
        for region_idx = 1:TARGET_REGIONS
            % Get target vertex ID for current brain region
            target_vertex_id = min_id(region_idx, 3);
            
            % Skip self-loop (source == target → shortest path = 0)
            if source_vertex_id ~= target_vertex_id
                % Compute shortest path between source and target vertex
                [~, shortest_path_length] = shortestpath(G, source_vertex_id, target_vertex_id);
                shortest_path_matrix(batch_vertex_count, region_idx) = shortest_path_length;
            else
                % Self-loop: set shortest path length to 0
                shortest_path_matrix(batch_vertex_count, region_idx) = 0;
            end
        end
    end
    
    %% ===================== Part 4: Save Batch Results =====================
    % Generate batch-specific filename (e.g., gd_1.mat, gd_2.mat)
    batch_str = num2str(batch_idx);
    result_filename = strcat('gd_', batch_str, '.mat');
    
    % Save results to dedicated directory (parallel-safe save)
    cd(RESULTS_DIR);
    parsave(result_filename, shortest_path_matrix);
    cd('..');  % Return to parent directory
    
    % Print batch completion status (progress tracking)
    fprintf('✅ Batch %d/%d completed! Results saved to: %s/%s\n', ...
        batch_idx, NUM_BATCHES, RESULTS_DIR, result_filename);
end

%% ===================== Part 5: Custom Parallel Save Function =====================
function parsave(fname, x)
    % parsave: Parallel-safe save function (avoids race conditions in batch saving)
    % Inputs:
    %   fname - Output filename (e.g., 'gd_1.mat')
    %   x - Shortest path matrix (60x91) for current batch
    save(fname, 'x');
end
