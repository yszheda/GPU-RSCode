bar(LogExp);
set(gca,'XTickLabel',{'approach 0', 'approach 1', 'approach 2', 'approach 3'});
legend('shared memory', 'constant memory', 'Location', 'Best');
ylabel('total GPU encoding time (ms)');