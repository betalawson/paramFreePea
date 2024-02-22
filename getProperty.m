function output = getProperty( input_array, property_name )
% This function extracts the property named <property_name> from each
% element in the array <input_array> and compiles them into a matrix, where
% each row is the property for that element. If the property being
% retrieved is matrix-valued, then this matrix is first reshaped into a 
% vector.

output = cellfun( @(x) x.(property_name)(:)', input_array(:), 'UniformOutput', false );

% Convert to a matrix if possible
try
    output = cell2mat(output);
end