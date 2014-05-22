figure;
plot(k2n3randomSpeedup(1,1:end), k2n3randomSpeedup(2,1:end));
set(gca,'XTick',1000:200:2000);
xlabel('file size (MB)');
ylabel('speedup');