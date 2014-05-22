plot(StepBW);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'GPU baseline implementation', 'optimize tiling algorithm', 'optimize log&exp table-based method', 'use pinned host memory', 'use CUDA streams'});
% legend('kernel execution', 'Location', 'Best');
ylabel('effective bandwidth (MB/s)');
rotateticklabel(gca,-30);

% get current (active) axes property
pos = get(gca,'OuterPosition');
bottomPadSize=0.2;
% set(gca,'OuterPosition',[left bottom + 0.1 width height]);
set(gca,'OuterPosition',[pos(1); pos(2) + bottomPadSize; pos(3) - 0.1; pos(4)-bottomPadSize]);