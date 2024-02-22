function K = kendallDist(X,Y)
% This function checks (in a slow, brute-force way) how many mismatched
% pairs exist in X and Y, when both are considered in an ordinal fashion.
% Ties are treated as wildcards (a tie between two elements in one list is
% satisfied regardless of the corresponding values in the other list)

% Read out length of lists
N = length(X);

% Initialise distance
K = 0;

if length(X) ~= length(Y)
    X
    Y
end

% Loop over all unique combinations of elements
for i = 1:N-1
    for j = i+1:N
        
        % Increment the distance counter in case of mismatch of rank
        if X(i) > X(j)
            if Y(i) < Y(j)
                K = K + 1;
            end
        elseif X(i) < X(j)
            if Y(i) > Y(j)
                K = K + 1;
            end
        end
        
    end
end

% Scale the value by the number of pairs (this makes lists of different
% lengths equally weighted)
K = K / (0.5*N*(N-1));