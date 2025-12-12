function print_metadata(structure, level)
    % Check if the input is a structure
    if ~isstruct(structure)
        disp('Input is not a structure');
        return;
    end

    % If no level is provided, start from level 0
    if nargin < 2
        level = 0;
    end
    
    % Get the field names of the structure
    fields = fieldnames(structure);
    
    % Loop through each field
    for i = 1:numel(fields)
        % Print indentation based on the level
        fprintf('%s%s: ', repmat('  ', 1, level), fields{i});
        
        % Check if the field contains another structure
        if isstruct(structure.(fields{i}))
            % If it's a structure, recursively call this function
%             fprintf('\n');
            print_metadata(structure.(fields{i}), level + 1);
        else
            % If it's not a structure, simply print its value
            disp(structure.(fields{i}));
        end
    end
end
