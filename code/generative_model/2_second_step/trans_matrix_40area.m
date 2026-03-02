% =========================================================================
% Function: trans_matrix_40area
% Purpose: Reduce a 91x91 neuron weight matrix to 40x40 by removing zero columns/rows
%          (zero columns are identified from MatF_SLN2021.mat)
% Inputs:
%   A       - Input matrix (typically 91x91 neuron weight matrix)
%   x       - Target number of rows (set to 40 to remove zero rows)
%   y       - Target number of columns (set to 40 to remove zero columns)
% Output:
%   B       - Reduced matrix (40x40 if x=40 and y=40; original A if x/y≠40)
% Dependencies:
%   - MatF_SLN2021.mat: Contains the reference 91x91 matrix to identify zero columns
% =========================================================================
function reduced_matrix = trans_matrix_40area(input_matrix, target_rows, target_cols)
    % Load reference matrix to identify zero columns (consistent with main analysis)
    load('MatF_SLN2021.mat');
    
    % Find columns with all zero values in the reference matrix
    zero_column_indices = find(all(MatF == 0, 1));
    
    % Create a copy of input matrix to avoid modifying original data
    reduced_matrix = input_matrix;
    
    % Remove zero rows (corresponding to zero columns) if target rows = 40
    if target_rows == 40 && ~isempty(zero_column_indices)
        reduced_matrix(zero_column_indices, :) = [];
    end
    
    % Remove zero columns if target columns = 40
    if target_cols == 40 && ~isempty(zero_column_indices)
        reduced_matrix(:, zero_column_indices) = [];
    end
    
   
end