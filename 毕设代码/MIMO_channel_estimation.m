clear; clc;close all; clf;
Nfft=512;
Ng=Nfft/8;
Nofdm=Nfft+Ng;
Nsym=10;
Nps=32; %导频间隔
Np=Nfft/Nps; %每个ofdm符号 导频数
Nd=Nfft-Np; % 每个ofdm符号 数据载波数
Nbps=4; % 每个映射符号的bit数
M=2^Nbps;
Es=1;% Signal energy
A=sqrt(3/2/(M-1)*Es); %  QAM normalization factor
SNRs = 20:20;  %dB
bers = zeros(3,length(SNRs));
MSE  = zeros(3,length(SNRs));
nose = zeros(3,length(SNRs));
Nt = 2;     % number of transmit antennas
Nr = 2;     % number of receive antenna

%%
for i=1:length(SNRs)
    SNR = SNRs(i);
    rng(12);  %rng 对应上层 每次仿真结果相同
    for nsym=1:Nsym
        rng(12);
        txSig=zeros(Nt,Nofdm);
        Xp=zeros(Nt,Np);
        XK=zeros(Nt,Nfft);
        pilot_loc=zeros(Nt,Np);
        msg=zeros(Nt,Nbps,Nfft-Np);
        for Nt_id=1:Nt
     
            Xp(Nt_id,:) = 2*(randn(1,Np)>0)-1;    % Pilot sequence generation
            msgint=round(rand(Nbps,Nfft-Np));     % bit generation
            msg(Nt_id,:,:)=msgint;
            Date = qammod(msgint,M,"gray","InputType","bit")*A;
           
            ip = 0;
            %pilot_loc =zeros(1,Np);
            X=zeros(1,Nfft);
            if Nt_id==1   %正交导频设计
                for k=1:Nfft
                    if mod(k,Nps)==1
                        X(k) = Xp(Nt_id,floor(k/Nps)+1);
                        pilot_loc(Nt_id,ip+1) = k;
                        ip = ip+1;
                    else
                        X(k) = Date(k-ip);
                    end
                    if mod(k,Nps)==2
                        X(k) = 0;
                    end
                end
                XK(Nt_id,:)=X;
                x = ifft(X,Nfft);                            % IFFT
                xt = [x(Nfft-Ng+1:Nfft) x];                  % Add CP
                txSig(Nt_id,:) =xt;        
            elseif Nt_id==2              
                for k=1:Nfft
                    if mod(k,Nps)==2
                        X(k) = Xp(Nt_id,floor(k/Nps)+1);
                        pilot_loc(Nt_id,ip+1) = k;
                        ip = ip+1;
                    else
                        X(k) = Date(k-ip);                       
                    end
                    if mod(k,Nps)==1
                        X(k) = 0;
                    end              
                end
                XK(Nt_id,:)=X;
                x = ifft(X,Nfft);                            % IFFT
                xt = [x(Nfft-Ng+1:Nfft) x];                  % Add CP
                txSig(Nt_id,:) =xt;           
            end         
        end
        %% 信道生成
        channel_length=10;
        h=zeros(Nr,Nt,channel_length);
        H=zeros(Nr,Nt,Nfft);
        pdp = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9 -25];% EVA
       %pdp=[0 -1.5];
        %power delay profile
        %pdp=zeros(1,10);
        pow_prof = 10.^(pdp/10);
        pow_prof = pow_prof/sum(pow_prof);%normalization of power delay profile
        rxSig=zeros(Nr,Nofdm+channel_length-1);
        for Nr_id=1:Nr           
            for Nt_id=1:Nt   % 抽头瑞利信道               
                chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1,channel_length)+1i*randn(1,channel_length)));%channel coef. for each path
          
                h(Nr_id,Nt_id,:)=chan_coef;
                temp_test1=squeeze(h(Nr_id,Nt_id,:));
                H(Nr_id,Nt_id,:)=fft(h(Nr_id,Nt_id,:),Nfft);
%                 temp_test2=squeeze(txSig(Nt_id,:));
%                                  hahahhah=temp_test2.';
                temp= conv(txSig(Nt_id,:),squeeze(h(Nr_id,Nt_id,:)));
%                                 H_test=toeplitz([temp_test1.'
%                                 zeros(1,574)],[temp_test1(1) zeros(1,574)
%                                 temp_test1(2)]); F=dftmtx(576); Fh=F';
%                                 F*H_test*Fh; y_test=H_test*temp_test2.';               
                %temp
                rxSig(Nr_id,:)=rxSig(Nr_id,:)+temp;
                rxSig(Nr_id,:)=awgn(rxSig(Nr_id,:),SNR,'measured');%前面是衰落，这里是加性噪声                
            end
        end
        rxSig=rxSig(:,Ng+1:Nofdm);
        rxSig=rxSig.';
        Y = fft(rxSig);
        Y=Y.';   % FFT  这里转换维度 因为fft的结果不同       
        H_LS=zeros(Nr,Nt,Nfft);
        H_MMSE=zeros(Nr,Nt,Nfft);
        H_LS_DFT=zeros(Nr,Nt,Nfft);
        H_MMSE_DFT=zeros(Nr,Nt,Nfft);
%         H_DFT_plus=zeros(Nr,Nt,Nfft);
%         H_DFT_plus_temp=zeros(1,Nfft);  
        for nr=1:Nr
            for nt=1:Nt         
                for m=1:2
                    if m==1             %信道估计
                        temp1=Y(nr,:); %因为不同维度的运算很麻烦，这里直接用temp存一维的数据
                        temp2=Xp(nt,:);
                        temp3=pilot_loc(nt,:);
                        temp4 = LS_CE(temp1,temp2,temp3,Nfft,Nps,'linear');                       
                        H_LS(nr,nt,:)=temp4;                        
                        h_est = ifft(temp4);                       
                        h_DFT = h_est(1:Ng);                               
                        H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation                       
                        H_LS_DFT(nr,nt,:)=H_DFT;   
                         %temp4_plus=LS_CE(temp1,temp2,temp3,Nfft,Nps,'spline');
                          %h_est_plus=ifft(temp4_plus);
                           % h_DFT_plus=h_est_plus(1:channel_length);            
%                         for index=1:size(H_DFT,2)
%                             if index==1
%                                 H_DFT_plus_temp(index)=H_DFT(1,index);
%                             elseif index>1
%                                 H_DFT_plus_temp(index)=H_DFT(1,index)/2+H_DFT(1,(index-1))/2;
%                             end
%                         end
                        %                                     H_DFT_plus_temp=fft(h_DFT_plus,Nfft);
                        %H_DFT_plus(nr,nt,:)=H_DFT_plus_temp;                      
                        temp7=squeeze(H(nr,nt,:)).'; %squeeze 去除1 维的维度                       
                        MSE(m,i) = MSE(m,i) + (temp7-temp4)*(temp7-temp4)';
                        MSE(m+1,i) = MSE(m+1,i) + (temp7-H_DFT)*(temp7-H_DFT)';                 
                        %MSE(4,i)=MSE(4,i)+(temp7-H_DFT_plus_temp)*(temp7-H_DFT_plus_temp)';
                    
                    elseif m==2
                        temp1=Y(nr,:);
                        temp2=Xp(nt,:);
                        temp3=pilot_loc(nt,:);
                        temp6=squeeze(h(nr,nt,:)).';
                        temp4 = MMSE_CE(temp1,temp2,temp3,Nfft,Nps,temp6,SNR);                        
                        H_MMSE(nr,nt,:)=temp4;                      
%                         h_est = ifft(temp4);
%                         h_DFT = h_est(1:channel_length);
%                         H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
%                         H_MMSE_DFT(nr,nt,:)=H_DFT;                        
                        %t=sum(abs(h_est(1:11)))/11+sum(abs(h_est(12:end)))/501;
                        %                                          t=sum(abs(h_est(1:11)))/100;
                        %                                          temp_DFT=h_est;
                        %                                          for
                        %                                          r=1:channel_length
                        %                                              if
                        %                                              abs(temp_DFT(r))<t
                        %                                                  temp_DFT(r)=0;
                        %                                              end
                        %
                        %                                          end
                        
                        %                                          H_DFT_plus(nr,nt,:)=fft(temp_DFT,Nfft);
                        %                                          H_DFT_plus_temp=
                        %                                          fft(temp_DFT,Nfft);
                                               
                        temp7=squeeze(H(nr,nt,:)).';                      
                        MSE(m+1,i) = MSE(m+1,i) + (temp7-temp4)*(temp7-temp4)';                       
                        %MSE(m+2,i) = MSE(m+2,i) +(temp7-H_DFT)*(temp7-H_DFT)';                       
                    end                  
                end               
            end
        end
        %%        
        %这里和SISO 信道估计不同 ，SISO估计只需要Y/H就行，MIMO需要进行联合运算才能解出X
        for m=1:2
            if m==1               
                H_tem=H_LS;
                X2 = (squeeze(H_tem(1,1,:)).'.*Y(2,:)-squeeze(H_tem(2,1,:)).'.*Y(1,:))./(squeeze(H_tem(1,1,:)).'.*squeeze(H_tem(2,2,:)).'-squeeze(H_tem(1,2,:)).'.*squeeze(H_tem(2,1,:)).');
                X1 =(Y(1,:)-squeeze(H_tem(1,2,:)).'.*X2)./squeeze(H_tem(1,1,:)).';
                ip = 0;
                Data_extracted=zeros(1,Nfft-Np);
                
                for k=1:Nfft
                    if mod(k,Nps)==1
                        ip=ip+1;
                    else
                        Data_extracted(k-ip)=X1(k);                       
                    end
                end
                figure(1)
%       plot(Data_extracted,'ro');
%       legend('rx symbs(LS)','rx symbs(MMSE)');
%     xlim([-2,2]); ylim([-2,2]);
%     grid on;
%     xlabel("Inphase");
%     ylabel("Quadrature");
                msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                temp_x1=squeeze(msg(1,:,:));
                nose(m,i) = nose(m,i) + sum(sum(msg_detected~=temp_x1));
                ip=0;
                Data_extracted=zeros(1,Nfft-Np);
                for k=1:Nfft
                    if mod(k,Nps)==2
                        ip=ip+1;
                    else
                        Data_extracted(k-ip)=X2(k);                       
                    end
                end
                msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                temp_x2=squeeze(msg(2,:,:));               
                nose(m,i) = nose(m,i) + sum(sum(msg_detected~=temp_x2));
                %%
                H_tem=H_LS_DFT;
                X2 = (squeeze(H_tem(1,1,:)).'.*Y(2,:)-squeeze(H_tem(2,1,:)).'.*Y(1,:))./(squeeze(H_tem(1,1,:)).'.*squeeze(H_tem(2,2,:)).'-squeeze(H_tem(1,2,:)).'.*squeeze(H_tem(2,1,:)).');
                X1 =(Y(1,:)-squeeze(H_tem(1,2,:)).'.*X2)./squeeze(H_tem(1,1,:)).';
                ip = 0;
                Data_extracted=zeros(1,Nfft-Np);
                for k=1:Nfft
                    if mod(k,Nps)==1
                        ip=ip+1;
                    else
                        Data_extracted(k-ip)=X1(k);                       
                    end
                end
                msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                temp_x1=squeeze(msg(1,:,:));
                nose(m+1,i) = nose(m+1,i) + sum(sum(msg_detected~=temp_x1));
                ip=0;
                Data_extracted=zeros(1,Nfft-Np);
                for k=1:Nfft
                    if mod(k,Nps)==2
                        ip=ip+1;
                    else
                        Data_extracted(k-ip)=X2(k);                        
                    end
                end
                msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                temp_x2=squeeze(msg(2,:,:));               
                nose(m+1,i) = nose(m+1,i) + sum(sum(msg_detected~=temp_x2));               
            elseif m==2
                H_tem=H_MMSE;
                X2 = (squeeze(H_tem(1,1,:)).'.*Y(2,:)-squeeze(H_tem(2,1,:)).'.*Y(1,:))./(squeeze(H_tem(1,1,:)).'.*squeeze(H_tem(2,2,:)).'-squeeze(H_tem(1,2,:)).'.*squeeze(H_tem(2,1,:)).');
                X1 =(Y(1,:)-squeeze(H_tem(1,2,:)).'.*X2)./squeeze(H_tem(1,1,:)).';
                ip = 0;
                Data_extracted=zeros(1,Nfft-Np);
                for k=1:Nfft
                    if mod(k,Nps)==1
                        ip=ip+1;
                    else
                        Data_extracted(k-ip)=X1(k);                       
                    end
                end
                 figure(1)
      plot(Data_extracted,'b^');
      legend('rx symbs(LS)','rx symbs(MMSE)');
    xlim([-2,2]); ylim([-2,2]);
    grid on;
    xlabel("Inphase");
    ylabel("Quadrature");
  
                msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                temp_x1=squeeze(msg(1,:,:));
                nose(m+1,i) = nose(m+1,i) + sum(sum(msg_detected~=temp_x1));
                ip=0;
                Data_extracted=zeros(1,Nfft-Np);
                for k=1:Nfft
                    if mod(k,Nps)==2
                        ip=ip+1;
                    else
                        Data_extracted(k-ip)=X2(k);                     
                    end
                end
                msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                temp_x2=squeeze(msg(2,:,:));                
                nose(m+1,i) = nose(m+1,i) + sum(sum(msg_detected~=temp_x2));
                %%
%                 H_tem=H_MMSE_DFT;
%                 X2 = (squeeze(H_tem(1,1,:)).'.*Y(2,:)-squeeze(H_tem(2,1,:)).'.*Y(1,:))./(squeeze(H_tem(1,1,:)).'.*squeeze(H_tem(2,2,:)).'-squeeze(H_tem(1,2,:)).'.*squeeze(H_tem(2,1,:)).');
%                 X1 =(Y(1,:)-squeeze(H_tem(1,2,:)).'.*X2)./squeeze(H_tem(1,1,:)).';
%                 ip = 0;
%                 Data_extracted=zeros(1,Nfft-Np);
%                 for k=1:Nfft
%                     if mod(k,Nps)==1
%                         ip=ip+1;
%                     else
%                         Data_extracted(k-ip)=X1(k);
%                         
%                     end
%                 end
                %msg_detected =
                %qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                %temp_x1=squeeze(msg(1,:,:)); nose(m+2,i) = nose(m+2,i) +
                %sum(sum(msg_detected~=temp_x1));
%                 ip=0;
%                 Data_extracted=zeros(1,Nfft-Np);
%                 for k=1:Nfft
%                     if mod(k,Nps)==2
%                         ip=ip+1;
%                     else
%                         Data_extracted(k-ip)=X2(k);
%                         
%                     end
%                 end
                %msg_detected =
                %qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
                %temp_x2=squeeze(msg(2,:,:));
                
                %nose(m+2,i) = nose(m+2,i) +
                %sum(sum(msg_detected~=temp_x2));
                    
            end

        end
        %真实的信道估计
%         H_tem=H;
%         X2 = (squeeze(H_tem(1,1,:)).'.*Y(2,:)-squeeze(H_tem(2,1,:)).'.*Y(1,:))./(squeeze(H_tem(1,1,:)).'.*squeeze(H_tem(2,2,:)).'-squeeze(H_tem(1,2,:)).'.*squeeze(H_tem(2,1,:)).');
%         X1 =(Y(1,:)-squeeze(H_tem(1,2,:)).'.*X2)./squeeze(H_tem(1,1,:)).';
%         ip = 0;
%         Data_extracted=zeros(1,Nfft-Np);
%         for k=1:Nfft
%             if mod(k,Nps)==1
%                 ip=ip+1;
%             else
%                 Data_extracted(k-ip)=X1(k);
%                 
%             end
%         end
%         msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
%         temp_x1=squeeze(msg(1,:,:));
        %nose(5,i) = nose(5,i) + sum(sum(msg_detected~=temp_x1));
        ip=0;
        Data_extracted=zeros(1,Nfft-Np);
        for k=1:Nfft
            if mod(k,Nps)==2
                ip=ip+1;
            else
                Data_extracted(k-ip)=X2(k);
                
            end
        end
%         msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
%         temp_x2=squeeze(msg(2,:,:));
        
        %nose(5,i) = nose(5,i) + sum(sum(msg_detected~=temp_x2));
        %改进的信道估计
%         H_tem=H_DFT_plus;
%         X2 = (squeeze(H_tem(1,1,:)).'.*Y(2,:)-squeeze(H_tem(2,1,:)).'.*Y(1,:))./(squeeze(H_tem(1,1,:)).'.*squeeze(H_tem(2,2,:)).'-squeeze(H_tem(1,2,:)).'.*squeeze(H_tem(2,1,:)).');
%         X1 =(Y(1,:)-squeeze(H_tem(1,2,:)).'.*X2)./squeeze(H_tem(1,1,:)).';
%         ip = 0;
%         Data_extracted=zeros(1,Nfft-Np);
%         for k=1:Nfft
%             if mod(k,Nps)==1
%                 ip=ip+1;
%             else
%                 Data_extracted(k-ip)=X1(k);
%                 
%             end
%         end
        %msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
%         temp_x1=squeeze(msg(1,:,:));
%         %nose(4,i) = nose(4,i) + sum(sum(msg_detected~=temp_x1));
%         ip=0;
%         Data_extracted=zeros(1,Nfft-Np);
%         for k=1:Nfft
%             if mod(k,Nps)==2
%                 ip=ip+1;
%             else
%                 Data_extracted(k-ip)=X2(k);
%                 
%             end
%         end
%        % msg_detected = qamdemod(Data_extracted/A,M,'gray','OutputType','bit');
%         temp_x2=squeeze(msg(2,:,:));
%         
        %nose(4,i) = nose(4,i) + sum(sum(msg_detected~=temp_x2));
     
    end
    
    MSE(:,i) = MSE(:,i)/(Nfft*Nsym*Nr*Nt);
    %    Number_of_symbol_errors= nose;
    for j=1:3
        bers(j,i)= nose(j,i)/((Nfft-Np)*Nsym*Nbps*Nt);
    end
%     
%     format compact %设置输出格式
%     clc;
%     display(bers(:,i),'BER');
%     display(MSE(:,i),'MSE');
%     display(SNRs(i),'SNR');
    
end
%% MSE
% figure(1)
% % LineStyles={'-o','-^','--o','--^','--s','--p'};
% % LineColor={[0 47 167]/255,[212 72 72]/255,[251 210 106]/255,[1 132 127]/255,[73 45 34]/255,[0 49 83]/255,[71 0 36]/255,[2 22 222]/255};
% LineStyles={'-o','-h','-^'};
% LineColor={[0 47 167]/255,[212 72 72]/255,'k'};
% for i=1:size(MSE,1)
%     semilogy(SNRs(1,:),MSE(i,:),LineStyles{i},'LineWidth',3,'MarkerSize',10,'Color',LineColor{i})
%     hold on
% end
% grid on
% xlabel('SNR/dB')
% ylabel('MSE')
% legend('LS','DFT转换域','MMSE')
% 
% title('2x2 MIMO-OFDM Channel Estimation,Nfft=512,16QAM')
% set(gca,'FontSize',20); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
% set(get(gca,'XLabel'),'FontSize',20);
% set(get(gca,'YLabel'),'FontSize',20);
% 
% %% BER
% figure(2)
% for i=1:size(bers,1)
%     semilogy(SNRs(1,:),bers(i,:),LineStyles{i},'LineWidth',3,'MarkerSize',10,'Color',LineColor{i})
%     hold on
% end
% grid on
% xlabel('SNR/dB')
% ylabel('BER')
% legend('LS','DFT转换域','MMSE')
% title('2x2 MIMO-OFDM Channel Estimation,Nfft=512,16QAM')
% set(gca,'FontSize',14); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
% set(get(gca,'XLabel'),'FontSize',20);%图上文字为8 point或小5号
% set(get(gca,'YLabel'),'FontSize',20);
% %%
% folderName = 'MIMO-OFDM'; % Folder name
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
