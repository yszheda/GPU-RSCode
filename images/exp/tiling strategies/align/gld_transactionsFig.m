figure;
bar(gld_transactions);
legend('original byte-length multiplication', 'word alignment','Location', 'Best');
xlabel('testcases');
ylabel('ratio of effective bandwidth compared to the highest one');