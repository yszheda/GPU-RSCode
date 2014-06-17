bar(LogExpNew);
set(gca,'XTickLabel',{'baseline', 'removing MOD(bit operations)', 'removing MOD(table augmentation)', 'further removing branches'});
% legend('shared memory (initialize serially)', 'constant memory', 'shared memory (load parallelly from constant memory)', 'Location', 'Best');
legend('constant memory', 'shared memory (initialize serially)', 'shared memory (load parallelly)', 'Location', 'Best');
ylabel('kernel execution time (ms)');
rotateticklabel(gca,-20);

% get current (active) axes property
pos = get(gca,'OuterPosition');
bottomPadSize=0.2;
% set(gca,'OuterPosition',[left bottom + 0.1 width height]);
set(gca,'OuterPosition',[pos(1); pos(2) + bottomPadSize; pos(3) - 0.1; pos(4) - bottomPadSize]);