function H_interpolated = interpolate(H_est,pilot_loc,Nfft,method,Nr,Nt)
% Input:        H_est    = Channel estimate using pilot sequence
%           pilot_loc    = location of pilot sequence
%                Nfft    = FFT size
%              method    = 'linear'/'spline'
% Output: H_interpolated = interpolated channel
%% 插值只会在已知区间内插值，所以要手动补全原点和终止点
% if pilot_loc(1)>1
%   slope = (H_est(2)-H_est(1))/(pilot_loc(2)-pilot_loc(1));
%   H_est = [H_est(1)-slope*(pilot_loc(1)-1)  H_est]; 
%   pilot_loc = [1 pilot_loc];
% end
% %%
% if pilot_loc(end)<Nfft
%     slope = (H_est(end)-H_est(end-1))/(pilot_loc(end)-pilot_loc(end-1));  
%     H_est = [H_est  H_est(end)+slope*(Nfft-pilot_loc(end))]; 
%     pilot_loc = [pilot_loc Nfft];
% end
% %%
% if lower(method(1))=='l'
%     H_interpolated = interp1(pilot_loc,H_est, 1:Nfft);   
% else
%     H_interpolated = interp1(pilot_loc,H_est, 1:Nfft ,'spline');
% end  

H_LS=zeros(Nr,Nt,Nfft);
for Nr_id=1:Nr
    for Nt_id=1:Nt
             %% 插值只会在已知区间内插值，所以要手动补全原点和终止点
             temp1=squeeze(pilot_loc(Nt_id,:));
             temp2=(squeeze(H_est(Nr_id,Nt_id,:))).';
            if pilot_loc(Nt_id,1)>1
              slope = (H_est(Nr_id,Nt_id,2)-H_est(Nr_id,Nt_id,1))/(pilot_loc(Nt_id,2)-pilot_loc(Nt_id,1));
              temp2 = [H_est(Nr_id,Nt_id,1)-slope*(pilot_loc(Nt_id,1)-1) temp2]; 
              temp1 = [1 temp1];
            end
            %%
            if pilot_loc(Nt_id,end)<Nfft
                slope = (H_est(Nr_id,Nt_id,end)-H_est(Nr_id,Nt_id,end-1))/(pilot_loc(Nt_id,end)-pilot_loc(Nt_id,end-1));  
                temp2 = [temp2  H_est(Nr_id,Nt_id,end)+slope*(Nfft-pilot_loc(Nt_id,end))]; 
               temp1 = [temp1 Nfft];
            end
            %%
            if lower(method(1))=='l'
                H_LS(Nr_id,Nt_id,:) = interp1(temp1,temp2, 1:Nfft);   
            else
                H_LS(Nr_id,Nt_id,:) = interp1(temp1,temp2, 1:Nfft ,'spline');
            end  


    end
end

H_interpolated=H_LS;