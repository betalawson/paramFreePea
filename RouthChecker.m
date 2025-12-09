function RouthChecker()

% This function calculates the Routh arrays using symbolic algebra
syms a b
syms m

%%% NO DELAY

% Routh array calculation - first two rows with zeros appended on end to ease below
R(1,:) = [1, 10, 5+2*a, 0, 0];
R(2,:) = [5, 10+a, 1+a+b, 0, 0];
% Routh array calculation - subsequent rows
for k = 3:6
    R(k,1:end-1) = ( R(k-1,1) .* R(k-2,2:end) - R(k-2,1) .* R(k-1,2:end) ) / R(k-1,1);
end

R_nodelay = simplify(R)



%%% DELAY

% Routh array calculation - first two rows with zeros appended on end to ease below
R(1,:) = [1, 21, 35+4*a, 7+4*a, 0];
R(2,:) = [7, 35+a, 21+6*a, 1+a+b, 0];
% Routh array calculation - subsequent rows
for k = 3:8
    R(k,1:end-1) = ( R(k-1,1) .* R(k-2,2:end) - R(k-2,1) .* R(k-1,2:end) ) / R(k-1,1);
end

R_delay = simplify(R)

% Consider the special case when b <= a
R_delay_no_fs = simplify(subs(R_delay, b, m*a))