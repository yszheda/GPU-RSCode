figure;
bar(k2n3randomStreamNum(1,1:end), k2n3randomStreamNum(2,1:end));
set(gca,'YTick',0:1:4);
xlabel('chunk size (MB)');
ylabel('best stream number');