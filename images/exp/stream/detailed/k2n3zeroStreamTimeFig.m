figure;
plot(k2n3zeroStreamTime);
legend('total GPU encoding time', 'kernel execution time', 'H2D data transfer time', 'Location', 'Best');
xlabel('number of CUDA streams');
ylabel('time (ms)');