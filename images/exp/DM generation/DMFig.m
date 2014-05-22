figure;
plot(DMTime(1,:), DMTime(2,:), 'Color', 'b', 'DisplayName', 'CPU');
hold on;
plot(DMTime(1,:), DMTime(3,:), 'Color', 'r', 'DisplayName', 'GPU');
hold off
legend('CPU', 'GPU', 'Location', 'Best');
xlabel('native chunk number (k)');
ylabel('time (ms)');