function [ output_args ] = p53Simulate( delta_t, T, initial )
%p53Simulate Runs simulation and outputs graph of p53 concentrations
%   Input initial concentrations?
%   initial - vector of initial concentrations.

results = zeros(T/delta_t+1, length(initial));
results(1, :) = initial;

time_vector = 0:delta_t:T;

for i = 2:T/delta_t+1
    
    %ODEs go here
    
end

%edit to only include p53 and modified p53 concentrations?
plot(time_vector, results);
title('p53 and Modified p53 Concentrations');

end

