%% MUSIC Doppler plots for each antenna and sum of all antennas
figure
for in = 1:12
    db2d = squeeze(db_psd_seg(in,:,:));
    subplot(3,4,in)
    pcolor(freq_seg,1:length(r),db2d)
    shading flat
    colormap('jet')
end
title([num2str(p_overlap_MU) '% overlap'])
drawnow

figure
db2d1 = movmean(db_psd_seg,fnpts,3);
db2d = squeeze(mean(db2d1,1));
pcolor(freq_seg,1:length(r),db2d)
shading flat
colormap('jet')
title(['antenna average ' num2str(fnpts) 'pt smooth with ' num2str(p_overlap_MU) '% overlap'])
drawnow