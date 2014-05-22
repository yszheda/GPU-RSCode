bar(timePerStep);
set(gca,'XTickLabel',{'baseline', 'optimize tiling algorithm', 'optimize log&exp table-based method', 'use pinned host memory', 'use CUDA streams'});
% legend('shared memory (initialize serially)', 'constant memory', 'shared memory (load parallelly from constant memory)', 'Location', 'Best');
% legend('shared memory (initialize serially)', 'constant memory', 'shared memory (load parallelly)', 'Location', 'Best');
ylabel('total GPU encoding time (ms)');
rotateticklabel(gca,30);

% get current (active) axes property
pos = get(gca,'OuterPosition');
bottomPadSize=0.2;
% set(gca,'OuterPosition',[left bottom + 0.1 width height]);
set(gca,'OuterPosition',[pos(1); pos(2) + bottomPadSize; pos(3); pos(4)-bottomPadSize]);