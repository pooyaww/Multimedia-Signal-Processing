%% MMSP2 - Lab 5
%  Exercise 1 - Audio coding with DPCM and ADPCM
%  Lucio Bianchi - 23/01/2014

clear
close all
clc

%% 1) Load the file '60.flac', convert to mono and take the segment from 1s to 3s
[x_stereo, Fs] = audioread('60.flac');
dur = 3 ;
x_left = x_stereo(:,1);
x_right = x_stereo(:,2);
x = (x_left+x_right)/2;
x = x(Fs:dur*Fs);

%% 2) PCM coding. Use an uniform quantizer with R = 2:16 bit. Estimate the SNR for each R.
%%    Plot the PSD of the input signal and of the reconstruction error
R = 2:16;
SNR_pcm = zeros(length(R),1);
x_max = 1;
x_min = -1;
for ii = 1:length(R)
    M = 2^(R(ii));
    Qx = (x_max - x_min)/M;
    x_tilde = Qx*floor(x/Qx);
    x_tilde(x < x_min) = x_min - Qx/2;
    x_tilde(x >= x_max) = x_max + Qx/2;
    SNR_pcm(ii) = db(var(x)/var(x-x_tilde));
    %sound(x_tilde,Fs)
%     disp(['Rate = ' num2str(R(ii))]);
%     nfft = 2^nextpow2(length(x));
%     Px = abs(fft(x,nfft)).^2;
%     Perr = abs(fft(x-x_tilde,nfft)).^2;
%     f_ax = linspace(0,Fs,nfft);
%     plot(f_ax(1:end/2),db(Px(1:end/2)),'b'), hold on
%     plot(f_ax(1:end/2),db(Perr(1:end/2)),'r'), hold off
%     legend('x(t)', 'error(t)')
%     title(['Power Spectral Density - Rate = ' num2str(R(ii))])
%     pause
end

%% 3) PCM coding and A-law companding. Repeat point a) and compare SNR. A = 87.7
SNR_a = zeros(length(R),1);
x_max = 1;
x_min = -1;
A = 87.7;
x_comp = zeros(length(x),1);
x_comp(abs(x) < 1/A) = A*x(abs(x)<1/A)/(1+log(A));
x_comp(abs(x) >= 1/A & abs(x) <= 1) = (1+log(A*abs(x(abs(x) >= 1/A & abs(x) <= 1))))/(1+log(A)).*sign(x(abs(x) >= 1/A & abs(x) <= 1));
for ii = 1:length(R)
    M = 2^R(ii);
    Qx = (x_max-x_min)/M;
    x_tilde = Qx*floor(x_comp/Qx);
    x_tilde(x_comp < x_min) = x_min - Qx/2;
    x_tilde(x_comp >= x_max) = x_max + Qx/2;
    x_exp = zeros(length(x_tilde),1);
    x_exp(abs(x_tilde) >= 0 & abs(x_tilde) <= 1/(1+log(A))) = x_tilde(abs(x_tilde) >= 0 & abs(x_tilde) <= 1/(1+log(A)))*(1+log(A))/A;
    x_exp(abs(x_tilde) > 1/(1+log(A)) & abs(x_tilde) <= 1) = exp(abs(x_tilde(abs(x_tilde) > 1/(1+log(A)) & abs(x_tilde) <= 1))*(1+log(A))-1)*1/A.*sign(x_tilde(abs(x_tilde) > 1/(1+log(A)) & abs(x_tilde) <= 1));
    SNR_a(ii) = db(var(x)/var(x-x_exp));
end

fig_snr = figure;
plot(R,SNR_pcm,'b'), hold on
plot(R,SNR_a,'r'), hold off
grid
xlabel('Rate'), ylabel('SNR [dB]')
legend('PCM','A-law','Location','BestOutside')

%% 4) DPCM. Repeat point a) and compare SNR. Predictor is x_hat(n) = alpha*x_tilde(n-1), alpha=0.9.
alpha = 0.9;
SNR_dpcm = zeros(length(R),1);
for ii = 1:length(R)
    M = 2^(R(ii));
    x_tilde = zeros(length(x),1);
    
    % first sample: PCM
    Qx = (x_max-x_min)/M;
    x_tilde(1) = Qx*floor(x(1)/Qx);
    
    % next samples: DPCM
    d_min = -0.1;
    d_max = 0.1;
    Qd = (d_max-d_min)/M;
    for nn = 2:length(x)
        d = x(nn) - alpha*x_tilde(nn-1);
        if d < d_min
            d_tilde = d_min - Qd/2;
        elseif d > d_max
            d_tilde = d_max + Qd/2;
        else
            d_tilde = Qd*floor(d/Qd);
        end
        x_tilde(nn) = d_tilde + alpha*x_tilde(nn-1);
    end
    SNR_dpcm(ii) = db(var(x)/var(x-x_tilde));
end

figure(fig_snr)
hold on
plot(R,SNR_dpcm,'g'), hold off
legend('PCM','A-law','DPCM','Location','BestOutside')

%% 5) APCM. Use Jayant multipliers for R = 2,3,4 and midrise quantizer.
R_apcm = 2:4;
SNR_apcm = zeros(length(R_apcm),1);
for ii = 1:length(R_apcm)
    M = 2^R_apcm(ii);
    switch M
        case 4
            Mult = [2.2 0.6 0.6 2.2];
        case 8
            Mult = [1.5 1 1 0.85 0.85 1 1 1.5];
        case 16
            Mult = [2.4 2 1.6 1.2 0.8 0.85 0.8 0.8 0.8 0.8 0.85 0.8 1.2 1.6 2 2.4];
    end
    
    Qx = zeros(length(x),1);
    x_tilde = zeros(length(x),1);
    % first sample: initialization
    Qx(1) = (x_max-x_min)/M;
    idx = -(M-1):2:(M-1);
    rec_levels = idx/2*Qx(1);
    x_tilde(1) = Qx(1)*floor(x(1)/Qx(1))+Qx(1)/2;
    [~,Mult_ind] = min(abs(x_tilde(1)-rec_levels));
    
    % next samples
    for nn = 2:length(x)
        Qx(nn) = Mult(Mult_ind) * Qx(nn-1);
        rec_levels = idx/2*Qx(nn);
        x_tilde(nn) = Qx(nn)*floor(x(nn)/Qx(nn))+Qx(nn)/2;
        [~,Mult_ind] = min(abs(x_tilde(nn)-rec_levels));
    end
    SNR_apcm(ii) = db(var(x)/var(x-x_tilde));
end

figure(fig_snr)
hold on
plot(R_apcm,SNR_apcm,'k'), hold off
legend('PCM','A-law','DPCM','APCM','Location','BestOutside')

%% 6) ADPCM. Repeat point a) and compare SNR. Predictor is x_hat(n) = alpha*x_tilde(n-1), alpha=0.9.
%%    Use Jayant multipliers for R = 2,3,4 and midrise quantizer.
alpha = 0.9;
R_adpcm = 2:4;
SNR_adpcm = zeros(length(R_adpcm),1);
for ii = 1:length(R_adpcm)
    M = 2^R_adpcm(ii);
    switch M
        case 4
            Mult = [1.6 0.8 0.8 1,6];
        case 8
            Mult = [1.75 1.25 0.9 0.9 0.9 0.9 1.25 1.75];
        case 16
            Mult = [2.4 2 1.6 1.2 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 1.2 1.6 2 2.4];
    end
    
    Qd = zeros(length(x),1);
    x_tilde = zeros(length(x),1);
    
    % first sample: PCM initialization
    Qx = (x_max-x_min)/M;
    idx = -(M-1):2:(M-1);
    rec_levels = idx/2*Qx;
    x_tilde(1) = Qx*floor(x(1)/Qx)+Qx/2;
    [~,Mult_ind] = min(abs(x_tilde(1)-rec_levels));
    Qd(1) = Qx;
    
    % next samples: DPCM
    for nn = 2:length(x)
        Qd(nn) = Mult(Mult_ind)*Qd(nn-1);
        rec_levels = idx/2*Qd(nn);
        d = x(nn) - alpha*x_tilde(nn-1);
        d_tilde = Qd(nn)*floor(d/Qd(nn));
        x_tilde(nn) = d_tilde + alpha*x_tilde(nn-1);
        [~,Mult_ind] = min(abs(d_tilde-rec_levels));
    end
    SNR_adpcm(ii) = db(var(x)/var(x-x_tilde));
end

figure(fig_snr)
hold on
plot(R_adpcm,SNR_adpcm,'m'), hold off
legend('PCM','A-law','DPCM','APCM','ADPCM','Location','BestOutside')