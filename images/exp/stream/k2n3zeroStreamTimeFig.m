figure;
plot(k2n3zeroStreamTime(1,:), k2n3zeroStreamTime(2,:), 'Color', 'b', 'DisplayName', 'total GPU encoding time (ms)');
% legend('total GPU encoding time (ms)', 'Location', 'Best');
xlabel('number of CUDA streams');
% ylabel('time (ms)');
ylabel('total GPU encoding time (ms)');