figure;
plot(k2n3randomStreamTime);
legend('total GPU encoding time', 'kernel execution time', 'non-overlapped data transfer time', 'Location', 'Best');
xlabel('number of CUDA streams');
ylabel('time (ms)');