% =========================================================================
% @function matrix_sort
% @brief   Custom matrix sorting function (sort by specified row/column)
% @author  [Ao Ma]
% @date    2026-03-02
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202531250002@mail.bnu.edu.cn]
% @affiliation [Beijing Normal University]
% 
% Description:
%   Sorts a matrix by a specified row or column in ascending/descending order
%   - Dim=1: Sort by COLUMN (row parameter = column index)
%   - Dim=2: Sort by ROW (row parameter = row index)
% 
% Examples:
%   1. Sort matrix A by 1st column (ascending): matrix_sort(A, 1, 1, 'ascend')
%   2. Sort matrix A by 2nd row (descending): matrix_sort(A, 2, 2, 'descend')
%   3. Default sort (1st column, ascending): matrix_sort(A, 1)
% =========================================================================

function [varargout] = matrix_sort(varargin)
    %% ===================== Input Validation & Defaults =====================
    % Check number of input parameters (2-4 required)
    if nargin < 2 || nargin > 4
        varargout{1} = {'Error: Missing or excess parameters! Required: 2-4 inputs'};
        return;
    end
    
    % Assign default values for optional parameters
    if nargin == 2
        varargin{3} = 1;          % Default: sort by column (dim=1)
        varargin{4} = 'ascend';   % Default: ascending order
    elseif nargin == 3
        varargin{4} = 'ascend';   % Default: ascending order
    end
    
    % Extract input parameters with clear naming
    input_matrix = varargin{1};
    target_index = varargin{2};
    sort_dimension = varargin{3};
    sort_order = varargin{4};
    
    % Validate sort dimension (must be 1 or 2)
    if sort_dimension ~= 1 && sort_dimension ~= 2
        varargout{1} = {'Error: Invalid 3rd parameter! Must be 1 (column) or 2 (row)'};
        return;
    end
    
    % Validate target index (within matrix bounds)
    if sort_dimension == 1  % Sort by column
        if target_index < 1 || target_index > size(input_matrix, 2)
            varargout{1} = {sprintf('Error: Column index %d out of bounds (matrix has %d columns)', target_index, size(input_matrix, 2))};
            return;
        end
    else  % Sort by row
        if target_index < 1 || target_index > size(input_matrix, 1)
            varargout{1} = {sprintf('Error: Row index %d out of bounds (matrix has %d rows)', target_index, size(input_matrix, 1))};
            return;
        end
    end
    
    % Validate sort order (must be 'ascend' or 'descend')
    if ~strcmp(sort_order, 'ascend') && ~strcmp(sort_order, 'descend')
        varargout{1} = {'Error: Invalid sort mode! Must be ''ascend'' or ''descend'''};
        return;
    end
    
    %% ===================== Matrix Sorting Logic =====================
    sorted_matrix = [];
    
    if sort_dimension == 1
        % Case 1: Sort by COLUMN (reorder rows based on target column)
        [~, sorted_indices] = sort(input_matrix(:, target_index), sort_order);
        sorted_matrix = input_matrix(sorted_indices, :);
    elseif sort_dimension == 2
        % Case 2: Sort by ROW (reorder columns based on target row)
        [~, sorted_indices] = sort(input_matrix(target_index, :), sort_order);
        sorted_matrix = input_matrix(:, sorted_indices);
    end
    
    %% ===================== Output Assignment =====================
    varargout{1} = sorted_matrix;
    
end