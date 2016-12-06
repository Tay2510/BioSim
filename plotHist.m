load('SSA.mat');

figure

before = p53helper_before./p53killer_before;
after = p53helper_after./p53killer_after;

[n1,x1] = hist(before);
h1 = bar(x1,n1/sum(n1),'hist'); hold on;
set(h1,'facecolor',[0 ,0.5, 0])

[n2,x2] = hist(after);
h2=bar(x2,n2/sum(n2),'hist');
set(h2,'facecolor','g')

box off
title('ratio distribution of p53helper/p53killer', 'fontsize', 16)
leg = legend('p53helper/p53killer','p53helper/p53killer (infected)');
set(leg,'fontsize', 16);
%xlabel('particle number', 'fontsize', 16)
xlabel('p53helper/p53killer', 'fontsize', 16)
ylabel('frequency', 'fontsize', 16)

disp('before')
disp(mean(before));
disp(var(before));
disp('after')
disp(mean(after));
disp(var(after));