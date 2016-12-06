num_trials = 1000;

p53_aggregate = zeros(num_trials, 2); %[post-virus p53_killer, post-virus p53_helper];

for i = 1:1000
    i
    p53_aggregate(i,:) = SDE();
end

edges = 0:.01:1;

histogram(p53_aggregate(:,1), edges);

figure 
histogram(p53_aggregate(:, 2), edges);