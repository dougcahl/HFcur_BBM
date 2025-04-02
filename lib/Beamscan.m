%% Beamscan
disp('calculating beamscan FFT...')
tic
[freq,zz,~,db_psd_ant,df] = fft_beam(I,Q,fs,w_BS,n_BS,p_overlap_BS);
freq = freq + f_offset; % add frequency offset
disp([num2str(toc) 's']) % time to run beamforming

%% Doppler plots for Beamform
if plot_doppler_beamscan == 1
    doppler_beamscan_plot;
end

%% 
disp('calculating Beamscan currents (DOAs) ...')
tic
% for each each range cel1, find the Bragg peaks, do MUSIC on them
for ir = 1:length(r)
    % get data for this range cell
    zphi1       = squeeze(a(:,:,ir));      % antenna pattern
    db1         =  squeeze(db_psd_ant(:,ir,:));  % fft db for each antenna, Beamscan

    % the matrices:  zz_mat , zz_seg_mat , zz_seg_coh_mat are passed directy to beamscan_range_ring
    zn = zz(ir,:); % antenna complex vals for Beamscan
    
    %% Beamscan
    [vdop00,vdop00n,theta00,theta00n,snr00,snr00n,signal00,signal00n] = ...
        beamscan_range_ring(zphi1,db1,zn,freq,f_bragg,lambda_bragg,vmax,db_min_BS,theta1,plot_db_peaks_BS,plot_BS_ind_spec,fnpts,bscan_usesmoothed_spec,w_bscan);

    beamscan_data(ir).vdop = vdop00;
    beamscan_data(ir).vdopn = vdop00n;
    beamscan_data(ir).theta = theta00;
    beamscan_data(ir).thetan = theta00n;
    beamscan_data(ir).snr = snr00;
    beamscan_data(ir).snrn = snr00n;
    beamscan_data(ir).signal = signal00;
    beamscan_data(ir).signaln = signal00n;
    
    %% plot Beamscan DOAs from this range cell
    if plot_beamscan_DOAs_per_range == 1
        sz = 50;
        figure
        subplot(311)
        [~,i_snr] = sort(snr00);
        scatter(theta00(i_snr),vdop00(i_snr),sz,snr00(i_snr),'filled','o')
        colorbar
        title(['r=' num2str(ir) ' positive peak - Beamscan'])

        subplot(312)
        [~,i_snr] = sort(snr00n);
        scatter(theta00n(i_snr),vdop00n(i_snr),sz,snr00n(i_snr),'filled','o')
        colorbar
        title('negative peak - Beamscan')

        subplot(313)
        [~,i_snr] = sort(snr00);
        scatter(theta00(i_snr),vdop00(i_snr),sz,snr00(i_snr),'filled','s')
        hold on
        [~,i_snr] = sort(snr00n);
        scatter(theta00n(i_snr),vdop00n(i_snr),sz,snr00n(i_snr),'filled','o')
        colorbar
        title('Beamscan - all DOAs')
        legend('pos peak','neg peak')
        drawnow

    end
end % beamscan on each range
disp([num2str(toc) 's'])

%% save data for this file
if save_beamscan == 1
    disp('saving Beamscan data...')
    save(fn_beamscan,'beamscan_data','theta_radar','r');
    disp('saved')
end

