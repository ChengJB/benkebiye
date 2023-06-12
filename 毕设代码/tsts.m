% MSE
SNR_CUT=SNRs(:,1:21);
MSE_CUT=MSE(:,1:21);
BSR_CUT=bers(:,1:21);
figure(1)
LineStyles={'-o','-h','-^'};
LineColor={[0 47 167]/255,[212 72 72]/255,'k'};
for i=1:size(MSE_CUT,1)
    if i<4
    semilogy(SNR_CUT(1,:),MSE_CUT(i,:),LineStyles{i},'LineWidth',3,'MarkerSize',10,'Color',LineColor{i})
    hold on
    end
end
grid on
xlabel('SNR/dB')
ylabel('MSE')
legend('LS','DFT转换域','MMSE')

title('2x2 MIMO-OFDM Channel Estimation,Nfft=512,16QAM')
set(gca,'FontSize',20); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);

%% BER
figure(2)
for i=1:size(BSR_CUT,1)
    
    if i<4
    semilogy(SNR_CUT(1,:),BSR_CUT(i,:),LineStyles{i},'LineWidth',3,'MarkerSize',10,'Color',LineColor{i})
    hold on
    end
end
grid on
xlabel('SNR/dB')
ylabel('BER')
legend('LS','DFT转换域','MMSE')
title('2x2 MIMO-OFDM Channel Estimation,Nfft=512,16QAM')
set(gca,'FontSize',14); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',20);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',20);

           