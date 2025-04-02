%% MUSIC
disp('calculating covariance matrices for MUSIC ...')
tic
% MUSIC covariance matrices
[freq_seg,db_psd_seg,zz_seg_mat,zz_seg_matM,zz_seg_matIM] = fft_cov_seg(I,Q,w_MU,fs,n_MU,p_overlap_MU);
freq_seg = freq_seg + f_offset; % add frequency offset
% hf_read_sorts_fft2; % freq_seg,db_psd_seg,zz_seg_mat,zz_seg_matM,zz_seg_matIM
disp([num2str(toc) 's'])
% clear I Q % to conserve ram

%% Doppler plots for MUSIC
if plot_doppler_MUSIC == 1
    doppler_MUSIC_plot;
end

%%
disp('calculating MUSIC currents (DOAs) ...')
tic
% for each each range cel1, find the Bragg peaks, do MUSIC on them
for ir = 1:length(r)
    % optional permuatation of 11 antennas or 10
    % optional comparison plots and calcs
    % plot range rings of radials currents, is it relatively smooth

    % get data for this range cell
    zphi1       = squeeze(a(:,:,ir));      % antenna pattern
    db1_MUSIC   =  squeeze(db_psd_seg(:,ir,:));  % fft db for antenna with overlapping segments, MUSIC

    % the matrices:  zz_mat , zz_seg_mat , zz_seg_coh_mat are passed directy to music_range_ring
    znn = []; % pass empty for MUSIC

    %% MUSIC
    % peak_thresh = 2;
    [vdop11,vdop11n,theta11,theta11n,snr11,snr11n,thetas11,thetas11n,signal11,signals11,signal11n,signals11n] = ...
        music_range_ring(zphi1,db1_MUSIC,znn,freq_seg,f_bragg,lambda_bragg,vmax,db_min_MUSIC,theta1,plot_db_peaks_MUSIC,plot_MUSIC_DOAs,zz_seg_mat,ir,fnpts,music_usesmoothed_spec,nmax,peak_thresh);
    music_data(ir).peak_thresh = peak_thresh;
    music_data(ir).vdop = vdop11;
    music_data(ir).vdopn = vdop11n;
    % music_data(ir).theta = theta11; % this is for a single peak (not accurate)
    % music_data(ir).thetan = theta11n;
    music_data(ir).thetas = thetas11;
    music_data(ir).thetasn = thetas11n;
    music_data(ir).snr = snr11;
    music_data(ir).snrn = snr11n;
    % music_data(ir).signal = signal11;
    % music_data(ir).signaln = signal11n;
    music_data(ir).signals = signals11;
    music_data(ir).signalsn = signals11n;

    %% Modified MUSIC for coherent signals, but the ocean is incoherent
    % MMUSIC
    % [vdop11,vdop11n,theta11,theta11n,snr11,snr11n,thetas11,thetas11n,signal11,signals11,signal11n,signals11n] = ...
    %     music_range_ring(zphi1,db1n,znn,freq_seg,f_bragg,lambda_bragg,vmax,db_min_MUSIC,theta1,n_dir,plot_db_peaks,plot_music00,vmax,zz_seg_matM,ir,fnpts,music_usesmoothed_spec,nmax,peak_thresh);
    % mmusic_data(ir).peak_thresh = peak_thresh;
    % mmusic_data(ir).vdop = vdop11;
    % mmusic_data(ir).vdopn = vdop11n;
    % mmusic_data(ir).theta = theta11;
    % mmusic_data(ir).thetan = theta11n;
    % mmusic_data(ir).thetas = thetas11;
    % mmusic_data(ir).thetasn = thetas11n;
    % mmusic_data(ir).snr = snr11;
    % mmusic_data(ir).snrn = snr11n;
    % mmusic_data(ir).signal = signal11;
    % mmusic_data(ir).signaln = signal11n;
    % mmusic_data(ir).signals = signals11;
    % mmusic_data(ir).signalsn = signals11n;

    % Improved Modified MUSIC for coherent signals, but the ocean is incoherent
    % IMMUSIC
    % [vdop11,vdop11n,theta11,theta11n,snr11,snr11n,thetas11,thetas11n,signal11,signals11,signal11n,signals11n] = ...
    %     music_range_ring(zphi1,db1n,znn,freq_seg,f_bragg,lambda_bragg,vmax,db_min_MUSIC,theta1,n_dir,plot_db_peaks,plot_music00,vmax,zz_seg_matIM,ir,fnpts,music_usesmoothed_spec,nmax,peak_thresh);
    % immusic_data(ir).peak_thresh = peak_thresh;
    % immusic_data(ir).vdop = vdop11;
    % immusic_data(ir).vdopn = vdop11n;
    % immusic_data(ir).theta = theta11;
    % immusic_data(ir).thetan = theta11n;
    % immusic_data(ir).thetas = thetas11;
    % immusic_data(ir).thetasn = thetas11n;
    % immusic_data(ir).snr = snr11;
    % immusic_data(ir).snrn = snr11n;
    % immusic_data(ir).signal = signal11;
    % immusic_data(ir).signaln = signal11n;
    % immusic_data(ir).signals = signals11;
    % immusic_data(ir).signalsn = signals11n;

    %% plot MUSIC DOAs from this range cell
    if plot_MUSIC_DOAs_per_range == 1
        sz = 50;
        figure
        subplot(311)
        [~,i_snr] = sort(signals11);
        if size(i_snr,1) > 1
            scatter(reshape(thetas11(i_snr),[],1),reshape(vdop11(i_snr),[],1),sz,reshape(signals11(i_snr),[],1),'filled','o')
        else
            scatter(reshape(thetas11(i_snr),[],1),repmat(vdop11,size(i_snr,2),1),sz,reshape(signals11(i_snr),[],1),'filled','o')
        end
        colorbar
        title('positive peak - MUSIC')

        subplot(312)
        [~,i_snr] = sort(signals11n);
        if size(i_snr,1) > 1
            scatter(reshape(thetas11n(i_snr),[],1),reshape(vdop11n(i_snr),[],1),sz,reshape(10*log10(signals11n(i_snr)),[],1),'filled','s')
        else
            scatter(reshape(thetas11n(i_snr),[],1),repmat(vdop11n,size(i_snr,2),1),sz,reshape(10*log10(signals11n(i_snr)),[],1),'filled','s')
        end
        colorbar
        title('negative peak - MUSIC')

        subplot(313)
        [~,i_snr] = sort(signals11);
        if size(i_snr,1) > 1
            scatter(reshape(thetas11(i_snr),[],1),reshape(vdop11(i_snr),[],1),sz,reshape(signals11(i_snr),[],1),'filled','o')
        else
            scatter(reshape(thetas11(i_snr),[],1),repmat(vdop11,size(i_snr,2),1),sz,reshape(signals11(i_snr),[],1),'filled','o')
        end
        hold on

        [~,i_snr] = sort(signals11n);
        if size(i_snr,1) > 1
            scatter(reshape(thetas11n(i_snr),[],1),reshape(vdop11n(i_snr),[],1),sz,reshape(10*log10(signals11n(i_snr)),[],1),'filled','s')
        else
            scatter(reshape(thetas11n(i_snr),[],1),repmat(vdop11n,size(i_snr,2),1),sz,reshape(10*log10(signals11n(i_snr)),[],1),'filled','s')
        end
        
        colorbar
        title('currents - MUSIC')
        legend('pos peak','neg peak')%,'beamform')
    end

    %% you can run MUSIC for 12 antennas, 11, 10, etc resulting in 100s of possible permutations
    % if do_permute == 1
    % end

end % music on each range
disp([num2str(toc) 's'])

%% save data for this file
if save_music == 1
    disp('saving MUSIC data...')
    save(fn_music,'music_data','theta_radar','r');
    disp('saved')
end
