function RouthCheckerDelayManual()

% Maximum values
a_max = 8;
b_max = 12;

% Loop points
Npts = 5001;

% Loop ranges
a_vals = linspace(0, a_max, Npts);
b_vals = linspace(0, b_max, Npts);

% Loop over these
for i = 1:Npts
    for j = 1:Npts
        
        a = a_vals(i);
        b = b_vals(j);
        
        % Only check scenarios that were unstable in the difference
        % equation model, these all involve either of root feedback signal
        % or root strigolactones turned off, which correspond to b <= a or
        % to b = 0, respectively, and the latter will also satisfy b <= a
        if b <= a
            
            % Routh array calculation - first two rows with zeros appended on end to ease below
            R(1,:) = [1, 21, 35+4*a, 7+4*a, 0];
            R(2,:) = [7, 35+a, 21+6*a, 1+a+b, 0];
            % Routh array calculation - subsequent rows
            for k = 3:8
                R(k,1:end-1) = ( R(k-1,1) .* R(k-2,2:end) - R(k-2,1) .* R(k-1,2:end) ) / R(k-1,1);
            end
            
            if ~all(R(:,1) > 0)
                keyboard
            end
        end
    end
end