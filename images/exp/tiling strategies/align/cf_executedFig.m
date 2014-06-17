figure;
bar(cf_executed);
legend('original byte-length multiplication', 'word alignment','Location', 'Best');
xlabel('testcases');
ylabel('ratio of effective bandwidth compared to the highest one');