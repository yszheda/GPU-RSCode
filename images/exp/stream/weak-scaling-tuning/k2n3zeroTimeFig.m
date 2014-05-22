figure;
plot(k2n3zeroTime(1,:), k2n3zeroTime(2,:), 'Color', 'b', 'DisplayName', 'H2D data transfer');
hold on;
plot(k2n3zeroTime(1,:), k2n3zeroTime(3,:), 'Color', 'r', 'DisplayName', 'encoding kernel');
hold off
legend('H2D data transfer', 'encoding kernel', 'Location', 'Best');
xlabel('file size(MB)');
ylabel('time (ms)');