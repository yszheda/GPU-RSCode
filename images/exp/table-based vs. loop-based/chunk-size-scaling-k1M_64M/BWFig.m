figure;
plot(chunkSizeScaling(1,2:end), chunkSizeScaling(2,2:end), 'Color', 'b', 'DisplayName', 'k=4,n=6(Table-based)');
hold on;
plot(chunkSizeScaling(1,2:end), chunkSizeScaling(3,2:end), 'Color', 'g', 'DisplayName', 'k=8,n=10(Table-based)');
plot(chunkSizeScaling(1,2:end), chunkSizeScaling(4,2:end), 'Color', 'm', 'DisplayName', 'k=16,n=18(Table-based)');
plot(chunkSizeScaling(1,2:end), chunkSizeScaling(5,2:end), 'Color', 'r', 'DisplayName', 'k=4,n=6(Loop-based)');
plot(chunkSizeScaling(1,2:end), chunkSizeScaling(6,2:end), 'Color', 'c', 'DisplayName', 'k=8,n=10(Loop-based)');
hold off
legend('k=4,n=6(Table-based)', 'k=8,n=10(Table-based)', 'k=16,n=18(Table-based)', 'k=4,n=6(Loop-based)', 'k=8,n=10(Loop-based)', 'Location', 'Best');
xlabel('chunk size (MB)');
ylabel('effective bandwidth (MB/s)');