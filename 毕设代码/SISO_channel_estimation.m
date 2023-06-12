%channel_estimation.m
% for LS/DFT Channel Estimation with linear/spline interpolation
clear; clc;close all; clf;
Nfft=512; 
Ng=Nfft/8;  
Nofdm=Nfft+Ng;
Nsym=1;
Nps=32; %导频间隔
Np=Nfft/Nps; %每个ofdm符号 导频数 
Nd=Nfft-Np; % 每个ofdm符号 数据载波数 
Nbps=4; % 每个映射符号的bit数
M=2^Nbps; 
Es=1;% Signal energy
A=sqrt(3/2/(M-1)*Es); %  QAM normalization factor
SNRs = 15:15;  %dB
bers = zeros(3,length(SNRs));
MSE =zeros(3,length(SNRs));
nose = zeros(3,length(SNRs));
% rng(1);
%%
for i=1:length(SNRs)
   SNR = SNRs(i); 
%    rng(12);
   for nsym=1:Nsym
      Xp = 2*(randn(1,Np)>0)-1;    % Pilot sequence generation
       
      msgint=round(rand(Nbps,Nfft-Np));    % bit generation
      Date = qammod(msgint,M,"gray","InputType","bit")*A;
      figure(1)
      %plot(Date,'r+');
      hold on;

      ip = 0;    
      pilot_loc =zeros(1,Np);
      X=zeros(1,Nfft);
      for k=1:Nfft
         if mod(k,Nps)==1
            X(k) = Xp(floor(k/Nps)+1);
            pilot_loc(1,ip+1) = k;
            ip = ip+1;
         else
             X(k) = Date(k-ip);
         end
      end
      x = ifft(X,Nfft);                            % IFFT
      xt = [x(Nfft-Ng+1:Nfft) x];                  % Add CP
      h = [(randn+1i*randn) (randn+1i*randn)/2 ] ;     % generates a (2-tap) channel
      H = fft(h,Nfft); 
      channel_length = length(h); % True channel and its time-domain length
      H_power_dB = 10*log10(abs(H.*conj(H)));      % True channel power in dB
      y_channel = conv(xt, h);                     % Channel path (convolution)
      %sig_pow = mean(y_channel.*conj(y_channel));
      yt = awgn(y_channel,SNR,'measured');  
      y = yt(Ng+1:Nofdm);                   %  remove delay spread (卷积长度会增加）and Remove CP
      Y = fft(y);                           %  FFT
      for m=1:2
         if m==1
             H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'linear'); 
             method='LS-linear'; % LS estimation with linear interpolation
             H_LS_power_dB = 10*log10(abs(H_est.*conj(H_est)));
%          elseif m==2
%              H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'spline');
%              method='LS-spline'; % LS estimation with spline interpolation
            H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
         h_est = ifft(H_est); 
         h_DFT = h_est(1:channel_length); 
         H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
        % H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));
         MSE(m+2,i) = MSE(m+2,i) + (H-H_DFT)*(H-H_DFT)';
         figure(1)
          
           plot(H_power_dB(1:15:end),'--o','linewidth',1);
            hold on;
           plot(H_LS_power_dB(1:15:end),'--^','linewidth',1);
           grid on; 
            xlabel('f/HZ')
            ylabel('power/dB')
legend('真实信道','LS估计信道')
title('MIMO-OFDM Channel Estimation ')
           
         else
             H_est = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR);
              H_MMSE_power_dB = 10*log10(abs(H_est.*conj(H_est)));
             
              figure(2)
               plot(H_power_dB(1:15:end),'--o','linewidth',1);
                hold on;
           plot(H_MMSE_power_dB(1:15:end),'--^','linewidth',1);
           xlabel('f/HZ')
            ylabel('power/dB')
legend('真实信道','MMSE估计信道')
title('MIMO-OFDM Channel Estimation ')
           %axis([0 32 -6 10])
           
           xlabel('Subcarrier Index'); 
           ylabel('Power[dB]');
             method='MMSE'; % MMSE estimation
         end
         
           
            
         MSE(m,i) = MSE(m,i) + (H-H_est)*(H-H_est)';
         
         %MSE(7,i) = MSE(7,i) + (H-H)*(H-H)';
 %%   没有DFT
      Y_eq = Y./H_est;
      ip = 0;
      Data_extracted=zeros(1,Nfft-Np);
      for k=1:Nfft
         if mod(k,Nps)==1
             ip=ip+1;  
         else
             Data_extracted(k-ip)=Y_eq(k); 
             
         end
      end
      msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
%       if m==1
%           figure(1)
%       plot(Data_extracted,'ro');
%       legend('rx symbs(LS)');
%     xlim([-2,2]); ylim([-2,2]);
%     grid on;
%     xlabel("Inphase");
%     ylabel("Quadrature");
%       elseif m==2
%           figure(1)
%       plot(Data_extracted,'b^');
%       legend('ref symbols','rx symbs(LSE)','rx symbs(MMSE)');
%    xlim([-2,2]); ylim([-2,2]);
%     grid on;
%     xlabel("Inphase");
%     ylabel("Quadrature");
%       end
      
      nose(m,i) = nose(m,i) + sum(sum(msg_detected~=msgint));
 %%  DFT
 if m==1
      Y_eq = Y./H_DFT;
      ip = 0;
      Data_extracted=zeros(1,Nfft-Np);
      for k=1:Nfft
         if mod(k,Nps)==1
             ip=ip+1;  
         else
             Data_extracted(k-ip)=Y_eq(k); 
             
         end
      end
      
      msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
      nose(m+2,i) = nose(m+2,i) + sum(sum(msg_detected~=msgint));
 end
      end
%%  理想信道
     Y_eq = Y./H;
      ip = 0;
      Data_extracted=zeros(1,Nfft-Np);
      for k=1:Nfft
         if mod(k,Nps)==1
             ip=ip+1;  
         else
             Data_extracted(k-ip)=Y_eq(k); 
             
         end
      end
      msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
%        figure(1)
%       plot(Data_extracted,'b^');
%       legend('rx symbs(MMSE)');
%    xlim([-2,2]); ylim([-2,2]);
%     grid on;
%     xlabel("Inphase");
%      ylabel("Quadrature");
%       nose(7,i) = nose(7,i) + sum(sum(msg_detected~=msgint));


   end   
   MSE(:,i) = MSE(:,i)/(Nfft*Nsym);
   Number_of_symbol_errors= nose;
   for j=1:size(bers,1)
       bers(j,i)= nose(j,i)/((Nfft-Np)*Nsym*Nbps);
   end

%    format compact %设置输出格式
%    clc;
%    display(bers(:,i),'BER');
%    display(MSE(:,i),'MSE');
%    display(SNRs(i),'SNR');

end   
%% MSE
% figure(1)
% LineStyles={'-*','-o','-^','-s','--o','--^','--s','-s'};
% LineColor={[0 47 167]/255,[212 72 72]/255,[251 210 106]/255,[1 132 127]/255,[73 45 34]/255,[0 49 83]/255,[71 0 36]/255,[2 22 222]/255};
% for i=1:size(MSE,1)
%     semilogy(SNRs(1,:),MSE(i,:),LineStyles{i},'LineWidth',2,'MarkerSize',8,'Color',LineColor{i})
%     hold on
% end
% grid on 
% xlabel('SNR/dB')
% ylabel('MSE')
% legend('LS','LS with DFT','MMSE')
% title('SISO-OFDM Channel Estimation 16QAM')
%% BER
% figure(2)
% for i=1:size(bers,1)
%     semilogy(SNRs(1,:),bers(i,:),LineStyles{mod(i,7)+1},'LineWidth',2,'MarkerSize',8,'Color',LineColor{mod(i,7)+1})
%     hold on
% end
% grid on 
% xlabel('SNR/dB')
% ylabel('BER')
% legend('LS','LS with DFT','MMSE')
% title('SISO-OFDM Channel Estimation 16QAM')
% %%
% folderName = 'test_1'; % Folder name
% month = datestr(now,5);date = datestr(now,7);
% folderPath=strcat('.\SimulationResult\',month,date,'\',folderName);
% if ~exist(folderPath,'dir')
%     mkdir(folderPath);%一运行就会在当前文件夹下创建simulation文件夹
% end 
% fileName = folderName; % File name
% %filePath = strcat(folderPath,'\',fileName,'.m');
% strsave = strcat(folderPath,'\',fileName,'.mat');
% s=['save ' strsave];% 保持.mat 文件，以后仿真结果可以再次确认,以后一定注意可以再次画图。
% eval(s);
