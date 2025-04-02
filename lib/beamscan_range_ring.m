function [vdop00,vdop00n,theta00,theta00n,snr00,snr00n,signal00,signal00n] = beamscan_range_ring(ant_pattern_phi,db_psd,zz,freq,f_bragg,lambda_bragg,vmax,db_min,theta1,plot_db_peaks,plot_BS_ind_spec,fnpts,bscan_usesmoothed_spec,w_beam)%#codegen
% antenner_pattern phi  =  is the phase response of the receiving antennas as a
% function of theta
% db_psd  =  fft db
% zz  =  fft complex coefficients
% freq  =  is the fft frequency vector
% f_bragg  =  is the bragg frequency
% lambda_bragg  =  is the Bragg wavelength
% vmax is the maximum radial speed in m/s
% db_min is the minimum db above noise to get music results
% theta1 is look direction from radar math notaion

% v_bragg_harmonic_diff
% for each frequency bin of interest, do music, save eigvals,
% snr, directions and p(theta) peak values and peaks
% fp0 = (freq - f_bragg); % Doppler from positive Bragg peak
% fn0 = (freq + f_bragg); % Dopper from negative Bragg peak

db_psd = nanmean(db_psd,1); % average antenna signal fft amplitudes
if bscan_usesmoothed_spec == 1
    db_psd = movmean(db_psd,fnpts,2);
end
db_psd = squeeze(mean(db_psd,1));
vp0 = (freq - f_bragg)*lambda_bragg; % velocities from positive Bragg peaks
vn0 = (freq + f_bragg)*lambda_bragg; % velocities from negative Bragg peaks
% v_bragg = f_bragg*lambda_bragg;
%%%%%%%%% for all frequencies within vmax, get average snr, do music db_min above this
fip0 = vp0 < vmax & vp0 > -vmax; % for positive Bragg peak in range +/- vmax
np0 = mean(db_psd(fip0)); % for positive Bragg peak(fip0)); % noise - data is above this

% find max and only within v_bragg_harmonic_diff and vmax
db_psda = db_psd;
db_psda(~fip0) = nan;
vpa = vp0;
vpa(~fip0) = nan;
[~,pmaxi] = max(db_psda);
% disp(['max dB p = ' num2str(pmax)])
vppeak = vp0(pmaxi);
% fip0 = abs(vppeak - vpa) < v_bragg_harmonic_diff/2;
% fip00 = fip0 & db_psd > np0 + db_min; % only peaks db_min above noise

fip00 = db_psda > np0 + db_min; % only peaks db_min above noise
freq00 = freq(fip00); % freq
vdop00 = vp0(fip00); % radial velocities
zz00 = zz(:,fip00); % complex z values
dbp00 = db_psda(fip00); % dB
snr00 = dbp00-np0; % snr above local Bragg peak noise noise

%%% neg freqs
fip0n = vn0 < vmax & vn0 > -vmax; % for positive Bragg peak in range +/- vmax
np0n = mean(db_psd(fip0n)); % for positive Bragg peak(fip0)); % noise - data is above this

% find max and only within v_bragg_harmonic_diff and vmax
db_psda = db_psd;
db_psda(~fip0n) = nan;
vna = vn0;
vna(~fip0n) = nan;
[~,pmaxi] = max(db_psda);
% disp(['max dB n = ' num2str(pmax)])
vnpeak = vn0(pmaxi);
% fip0n = abs(vnpeak - vna) < v_bragg_harmonic_diff/3;
% fip00n = fip0n & db_psd > np0n + db_min; % only peaks db_min above noise

fip00n = db_psda > np0n + db_min; % only peaks db_min above noise
freq00n = freq(fip00n); % freq
vdop00n = vn0(fip00n); % radial velocities
zz00n = zz(:,fip00n); % complex z values
dbp00n = db_psda(fip00n); % dB
snr00n = dbp00n-np0n; % snr above noise

% optional plot
if plot_db_peaks == 1
    figure
    plot(freq,db_psd)
    hold on
    plot(freq00,dbp00,'o')
    plot(freq00n,dbp00n,'o')
    a = ylim;
    plot([f_bragg f_bragg],a,':k')
    plot([-f_bragg -f_bragg],a,':k')
    title('Beamscan frequencies used for DOAs')
    drawnow
%     pause
end

%% Scan over theta
[theta00,signal00] = beamscan_single_bragg_peakc(ant_pattern_phi,freq00,zz00,theta1,plot_BS_ind_spec,'positive peak',w_beam);
[theta00n,signal00n] = beamscan_single_bragg_peakc(ant_pattern_phi,freq00n,zz00n,theta1,plot_BS_ind_spec,'negative peak',w_beam);

% mex files don't do multiple DOAs correctly? They are different than the m
% file above
% [atheta00,athetas00,asignal00,asignals00] = music_single_bragg_peakc_mex(ant_pattern_phi,freq00,zz00,theta1,n_dir,plot_music00,'positive peak',zzmat2);
% [theta00n,thetas00n,signal00n,signals00n] = music_single_bragg_peakc_mex(ant_pattern_phi,freq00n,zz00n,theta1,n_dir,plot_music00,'negative peak',zzmat2n);
