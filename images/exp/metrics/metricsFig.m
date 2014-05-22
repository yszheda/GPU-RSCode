% plot(metrics(:,1), metrics(:,2), '+');
% xlabel('kernel execution time (ms)');
% title('cf_executed');
% saveas(gcf,'cf_executed.pdf');
% close all;
% 
% plot(metrics(:,1), metrics(:,3), '+');
% xlabel('kernel execution time (ms)');
% title('cf issued');
% saveas(gcf,'cf_issued.pdf');
% close all;
% 
% plot(metrics(:,1), metrics(:,4), '+');
% xlabel('kernel execution time (ms)');
% title('dram read transactions');
% saveas(gcf,'dram_read_transactions.pdf');
% close all;
% 
% plot(metrics(:,1), metrics(:,5), '+');
% xlabel('kernel execution time (ms)');
% title('gld efficiency');
% saveas(gcf,'gld_efficiency.pdf');
% close all;
% 
% plot(metrics(:,1), metrics(:,6), '+');
% xlabel('kernel execution time (ms)');
% title('gld transactions');
% saveas(gcf,'gld_transactions.pdf');
% close all;

for i = 2:35,
    plot(metrics(:,1), metrics(:,i), '+');
    xlabel('kernel execution time (ms)');
    title(metrics_name(i));
    fileName = char(strcat(metrics_name(i), '.pdf'));
    saveas(gcf, fileName);
    close all;
end
