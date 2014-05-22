figure;
plot(k2n3zeroStreamSpeedup(1,:), k2n3zeroStreamSpeedup(end,:), 'Color', 'b', 'DisplayName', 'speedup');
% plot(k2n3zeroStreamSpeedup(end,:), 'Color', 'b', 'DisplayName', 'speedup');
% plot(k2n3zeroStreamSpeedup(end,:));
% legend('speedup', 'Location', 'Best');
xlabel('number of CUDA streams');
% ylabel('time (ms)');
ylabel('speedup');