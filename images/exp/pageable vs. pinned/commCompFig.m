bar(memComp);
set(gca,'XTickLabel',{'H2D Data Transfer', 'D2H Data Transfer'});
legend('pageable memory', 'pinned memory', 'Location', 'Best');
ylabel('time (ms)');