%% Beamform
%% first calculate FFT then beamform
disp('calculating beamforming FFT...')
tic
% freq, zz (covariance matrices R), db_psd
[freq,z_beamform,db_psd,db_psd_ant,df] = fft_beam(I,Q,fs,w_BF,n_BF,p_overlap_BF);
freq = freq + f_offset; % add frequency offset
disp([num2str(toc) 's']) % time to run beamforming

%% Beam straight forward pcolor
if plot_doppler_Bform == 1
    figure
    pcolor(freq,r,db_psd)
    colormap(jet)
    shading flat
    colorbar
    title('all antennas summed')
    xlabel('Freq (Hz)')
    ylabel('Range (km)')
    drawnow
end

%% Plot all antennas
if plot_doppler_Bform_ants == 1
    for iplt = 1:antN
        figure
        pcolor(freq,r/10^3,squeeze(db_psd_ant(iplt,:,:)))
        colormap(jet)
        shading flat
        colorbar
        title(['antenna #' num2str(iplt)])
        xlabel('Freq (Hz)')
        ylabel('Range (km)')
        drawnow
    end
end

%% Beamforming
disp('calculating Doppler beamformed spectra ...')
tic
zbeam = beamform11c(z_beamform,a,w_bform);
db_beam = 10*log10(abs(zbeam)/df);
disp([num2str(toc) 's']) % time to run beamforming

%% Beam straight forward pcolor
if plot_doppler_Bform == 1
    figure
    pcolor(freq,r/10^3,squeeze(db_beam(:,:,ceil(size(db_beam,3)/2))))
    shading flat
    colormap(jet)
    colorbar
    title(['Beamformed to ' num2str(theta1(ceil(size(db_beam,3)/2))) '^\circ'])
    xlabel('Freq (Hz)')
    ylabel('Range (km)')
    drawnow

    % Beam left
    figure
    pcolor(freq,r/10^3,squeeze(db_beam(:,:,1)))
    shading flat
    colormap(jet)
    colorbar
    title(['Beamformed to ' num2str(theta1(1)) '^\circ'])
    xlabel('Freq (Hz)')
    ylabel('Range (km)')
    drawnow

    % Beam right
    figure
    pcolor(freq,r/10^3,squeeze(db_beam(:,:,size(zbeam,3))))
    shading flat
    colormap(jet)
    colorbar
    title(['Beamformed to ' num2str(theta1(end)) '^\circ'])
    xlabel('Freq (Hz)')
    ylabel('Range (km)')
    drawnow
end

%% calculate current from biggest Bragg peak and Doppler spec of
% beamformed SORT data
disp('calculating surface current estimates ...')
tic

V = nan(size(zbeam,1),size(zbeam,3));
SNR = V;
NOISE = V;
itheta1 = 1;
for itheta = 1:dtheta_beam:length(theta1)
    db_psd1 = squeeze(db_beam(:,:,itheta));

    % using the mex file will go faster but you should compile
    % your own mex file
    [v,snr,noise_beam] = hf_beam_max_current_upd1(freq,db_psd1,vmax,lambda_bragg,f_bragg,npts,plot_beamform_cur,theta1(itheta));

    V(:,itheta1) = v;
    SNR(:,itheta1) = snr;
    NOISE(:,itheta1) = noise_beam;
    itheta1 = itheta1 + 1;
end
disp([num2str(toc) 's'])

%% save data
if save_data_beamform == 1 % FFT spectra and beamformed spectra
    save(fn_beamIQ,'theta1','z_beamform','zbeam','theta_radar','r')%,'theta_beam')
end

if save_data_beamform_cur == 1 % beamformed currents
    disp('saving data ...')
    tic
    save(fn_beamform,'V','SNR','NOISE','theta1','theta_radar','r')%,'theta_beam')
    disp([num2str(toc) 's'])
end

