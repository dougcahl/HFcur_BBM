%% Beamscan Doppler plots for each antenna and sum of all antennas
Nbeamscan = size(db_psd_ant,3);

figure
for in = 1:12
    db2d = squeeze(db_psd_ant(in,:,:));
    subplot(3,4,in)
    pcolor(freq,1:length(r),db2d)
    shading flat
    colormap('jet')
end
title(['full ' num2str(Nbeamscan) ' pt FFT'])
drawnow

figure
for in = 1:12
    db2d = squeeze(db_psd_ant(in,:,:));
    subplot(3,4,in)
    pcolor(freq,1:length(r),movmean(db2d,fnpts,2))
    shading flat
    colormap('jet')
end
title([num2str(fnpts) 'pt smooth full ' num2str(Nbeamscan) 'pt FFT'])
drawnow

