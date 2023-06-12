function H_LS = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,int_opt,Nr,Nt)
% LS channel estimation function
% Inputs:
%       Y         = Frequency-domain received signal
%       Xp        = Pilot signal
%       pilot_loc = Pilot location
%       N         = FFT size
%       Nps       = Pilot spacing
%       int_opt   = 'linear' or 'spline'
% output:
%       H_LS      = LS channel etimate
Np=Nfft/Nps;
LS_est=zeros(Nr,Nt,Np);
H_LS=zeros(Nr,Nt,Nfft);

for i=1:Nr
    for j=1:Nt
        for k=1:Np
            LS_est(i,j,k) = Y(i,pilot_loc(j,k))./Xp(j,k);  % LS channel estimation
            
        end


    end
end

if  lower(int_opt(1))=='l'
    method='linear'; 
else
    method='spline'; 
end

H_LS = interpolate(LS_est,pilot_loc,Nfft,method,Nr,Nt); % Linear/Spline interpolation
% for k=1:Np
%  LS_est(k) = Y(pilot_loc(k))./Xp(k);  % LS channel estimation
% end

% if  lower(int_opt(1))=='l'
%     method='linear'; 
% else
%     method='spline'; 
% end
% H_LS = interpolate(LS_est,pilot_loc,Nfft,method); % Linear/Spline interpolation