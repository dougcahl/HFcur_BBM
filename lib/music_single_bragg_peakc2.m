function [theta00,thetas00,signal1,signals1] = music_single_bragg_peakc2(ant_pattern_phi,freq00,theta1,plot_MUSIC_DOAs,tname,zzmat2,nmax,peak_thresh)

theta00 = nan(length(freq00),1);
thetas00 = nan(length(freq00),nmax);
signal1 = nan(length(freq00),1);
signals1 = nan(length(freq00),nmax);
if plot_MUSIC_DOAs == 1
    figure
end
for i_00 = 1:length(freq00)
    if plot_MUSIC_DOAs == 1
        clf
        title(tname)
    end
    Z = zzmat2{i_00};
    [atheta00,athetas00,asignal1,asignals1] = music_single_bragg_peakc1_sub(ant_pattern_phi,theta1,nmax,peak_thresh,Z,plot_MUSIC_DOAs,tname);
    
    thetas00(i_00,:) = athetas00;
    signals1(i_00,:) = asignals1; 
    
    theta00(i_00) = atheta00;
    signal1(i_00) = asignal1; 

end % for each frequency peak

