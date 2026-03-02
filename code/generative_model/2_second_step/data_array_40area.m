% =========================================================================
% @function data_array_40area
% @brief   Extract and filter 40-area feature array from 91x40/91x91 neuron matrix
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% 
% Description:
%   1. Flattens input matrix (A) into a 3-column feature array (row_idx, col_idx, value)
%   2. Filters out zero columns from MatF_SLN2021.mat to get 40 valid brain areas
%   3. Removes specific index pairs (ind91(1,i), ind91(2,i)) from the feature array
% 
% =========================================================================

function Data_A = data_array_40area(input_matrix)
    %% ===================== Input Validation =====================
    % Check input matrix dimensions (must be 91x40 or 91x91 for neuron analysis)
    matrix_size = size(input_matrix);
    if ~( (matrix_size(1) == 91 && matrix_size(2) == 40) || (matrix_size(1) == 91 && matrix_size(2) == 91) )
        error('data_array_40area:InvalidMatrixSize', ...
              'Input matrix must be 91x40 or 91x91 (received %dx%d)', matrix_size(1), matrix_size(2));
    end
    
    %% ===================== Step 1: Flatten Input Matrix to Feature Array =====================
    % Initialize feature array (3 columns: row index, column index, matrix value)
    total_elements = matrix_size(1) * matrix_size(2);
    Data_A = zeros(total_elements, 3);
    element_count = 0;
    
    % Iterate through matrix and populate feature array
    for row_idx = 1:matrix_size(1)
        for col_idx = 1:matrix_size(2)
            element_count = element_count + 1;
            Data_A(element_count, 1) = row_idx;          % Row index in input matrix
            Data_A(element_count, 2) = col_idx;          % Column index in input matrix
            Data_A(element_count, 3) = input_matrix(row_idx, col_idx);  % Matrix value (gradient/weight)
        end
    end
    
    %% ===================== Step 2: Load Weight Matrix & Filter Valid 40 Areas =====================
    % Load MatF_SLN2021.mat to identify zero columns (invalid brain areas)
    load('MatF_SLN2021.mat');
    
    % Find zero columns (columns with all zero values in MatF)
    zero_column_indices = find(all(MatF == 0, 1));
    
    % Create index list for 91 areas, remove zero columns (retains 40 valid areas)
    full_area_indices = 1:91;
    valid_area_indices_91 = full_area_indices;
    valid_area_indices_91(zero_column_indices) = [];
    
    % Create 2-row index matrix (row1=valid 91-area indices, row2=1:40 mapped indices)
    valid_area_mapping = zeros(2, 40);
    valid_area_mapping(1, :) = valid_area_indices_91;  % Row 1: valid 91-area indices
    valid_area_mapping(2, :) = 1:40;                   % Row 2: 1-40 mapped indices
    
    %% ===================== Step 3: Remove Specific Index Pairs =====================
    % Find indices to delete (matching (ind91(1,i), ind91(2,i)) pairs)
    indices_to_delete = [];
    for area_idx = 1:40
        % Get target index pair for current valid area
        target_row = valid_area_mapping(1, area_idx);
        target_col = valid_area_mapping(2, area_idx);
        
        % Find matching rows in Data_A (vectorized for speed)
        match_mask = (Data_A(:,1) == target_row) & (Data_A(:,2) == target_col);
        match_indices = find(match_mask);
        
        % Accumulate indices to delete (avoid empty/duplicate indices)
        if ~isempty(match_indices)
            indices_to_delete = [indices_to_delete; match_indices];
        end
    end
    
    % Remove duplicate indices (edge case: multiple matches)
    indices_to_delete = unique(indices_to_delete);
    
    % Delete specified rows from Data_A (only if indices exist)
    if ~isempty(indices_to_delete)
        Data_A(indices_to_delete, :) = [];
    else
        warning('data_array_40area:NoIndicesToDelete', ...
                'No matching index pairs found for deletion (all 40-area entries retained)');
    end
    
    %% ===================== Final Sanity Check =====================
    % Remove any all-zero rows (edge case: empty values after filtering)
    Data_A = Data_A(~all(Data_A == 0, 2), :);
    
end
