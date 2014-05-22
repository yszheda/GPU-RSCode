figure;
 plot(k2n3randomStreamScaleSpeedup(:,1), k2n3randomStreamScaleSpeedup(:,3:end));
legend('#stream=2', '#stream=3', '#stream=4', 'Location', 'Best');
set(gca,'XTick',1000:200:2000);
xlabel('file size (MB)');
ylabel('speedup');