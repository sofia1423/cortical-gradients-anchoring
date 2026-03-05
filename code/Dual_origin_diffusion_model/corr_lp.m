% =========================================================================
% @file    corr_lp.m
% @brief   Calculate Pearson Correlation for Non-Zero Pairs of Two Vectors
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% @note    Key logic: Filter (lp, lp1) pairs where both values ≠ 0 → compute Pearson correlation
% =========================================================================

function [corr_coeff, p_value] = corr_lp(input_vector1, input_vector2)
    %% ===================== Input Validation =====================
    % Check input dimensions (must be column vectors with same length)
    if size(input_vector1, 2) ~= 1 || size(input_vector2, 2) ~= 1
        error('corr_lp:InvalidInputDimension', ...
              'Inputs must be column vectors (received input1: %dx%d, input2: %dx%d)', ...
              size(input_vector1,1), size(input_vector1,2), size(input_vector2,1), size(input_vector2,2));
    end
    
    if length(input_vector1) ~= length(input_vector2)
        error('corr_lp:InputLengthMismatch', ...
              'Input vectors must have same length (input1: %d, input2: %d)', ...
              length(input_vector1), length(input_vector2));
    end
    
    %% ===================== Filter Non-Zero Pairs =====================
    % Create logical mask for non-zero pairs (both input1 and input2 ≠ 0)
    non_zero_mask = (input_vector1 ~= 0) & (input_vector2 ~= 0);
    
    % Count valid non-zero pairs
    num_valid_pairs = sum(non_zero_mask);
    
    % Handle edge case: no valid non-zero pairs
    if num_valid_pairs < 2  % Need at least 2 pairs for correlation
        warning('corr_lp:InsufficientValidPairs', ...
                'Only %d valid non-zero pairs (need ≥2) → returning corr=NaN, p=NaN', num_valid_pairs);
        corr_coeff = NaN;
        p_value = NaN;
        return;
    end
    
    % Extract valid non-zero pairs (vectorized - faster than loop)
    valid_vector1 = input_vector1(non_zero_mask);
    valid_vector2 = input_vector2(non_zero_mask);
    
    %% ===================== Calculate Pearson Correlation =====================
    % Compute Pearson correlation (same as original corr() call)
    [corr_coeff, p_value] = corr(valid_vector1, valid_vector2);
    
    %% ===================== Optional Debug Info =====================
    % Uncomment below to print validation stats (helpful for debugging)
    % fprintf('corr_lp completed: %d valid non-zero pairs | Correlation: %.4f | p-value: %.4f\n', ...
    %     num_valid_pairs, corr_coeff, p_value);
end