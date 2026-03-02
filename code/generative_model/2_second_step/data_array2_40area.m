% =========================================================================
% @function data_array2_40area
% @brief   Extract non-zero 40-area feature array from neuron gradient/weight matrix
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% 
% Description:
%   1. Transforms MatF (91x91 weight matrix) to 91x40 matrix using trans_matrix_40area
%   2. Extracts non-zero entries from input matrix (A) where M91n (91x40 MatF) is non-zero
%   3. Returns 3-column feature array (row_idx, col_idx, value) for valid non-zero entries
% =========================================================================

function Data_A = data_array2_40area(input_matrix)
    %% ===================== Input Validation =====================
    % Check input matrix dimensions (must be 91x40 or 91x91 for neuron analysis)
    matrix_size = size(input_matrix);
    if ~( (matrix_size(1) == 91 && matrix_size(2) == 40) || (matrix_size(1) == 91 && matrix_size(2) == 91) )
        error('data_array2_40area:InvalidMatrixSize', ...
              'Input matrix must be 91x40 or 91x91 (received %dx%d)', matrix_size(1), matrix_size(2));
    end
    
    %% ===================== Step 1: Load & Transform Weight Matrix to 91x40 =====================
    % Load core weight matrix (MatF_SLN2021.mat)
    load('MatF_SLN2021.mat');
    
    % Transform 91x91 MatF to 91x40 valid area matrix (M91n)
    weight_matrix_91x40 = trans_matrix_40area(MatF, 91, 40);
    
    %% ===================== Step 2: Extract Non-Zero Entries (M91n != 0) =====================
    % Precompute non-zero mask (vectorized for speed)
    non_zero_mask = weight_matrix_91x40 ~= 0;
    total_non_zero = sum(non_zero_mask(:));
    
    % Validate non-zero entries (avoid empty output)
    if total_non_zero == 0
        warning('data_array2_40area:NoNonZeroEntries', ...
                'No non-zero entries found in 91x40 weight matrix (weight_matrix_91x40)');
        Data_A = [];  % Return empty array if no valid entries
        return;
    end
    
    % Pre-allocate feature array (faster than dynamic growth)
    Data_A = zeros(total_non_zero, 3);
    element_count = 0;
    
    % Iterate through matrix and extract non-zero entries (M91n != 0)
    for row_idx = 1:matrix_size(1)
        for col_idx = 1:matrix_size(2)
            % Only include entries where weight_matrix_91x40 is non-zero
            if non_zero_mask(row_idx, col_idx)
                element_count = element_count + 1;
                Data_A(element_count, 1) = row_idx;          % Row index in input matrix
                Data_A(element_count, 2) = col_idx;          % Column index in input matrix
                Data_A(element_count, 3) = input_matrix(row_idx, col_idx);  % Input matrix value
            end
        end
    end
    
    %% ===================== Final Sanity Check =====================
    % Ensure pre-allocation matches actual count (edge case handling)
    if element_count < total_non_zero
        Data_A = Data_A(1:element_count, :);  % Trim excess pre-allocated rows
    end
    
end
