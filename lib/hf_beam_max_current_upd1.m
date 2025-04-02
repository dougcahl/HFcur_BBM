function [v,snr,noise_beam] = hf_beam_max_current_upd1(freq,db_psd,vmax,lambda_bragg,f_bragg,npts,plot_beamform_cur,theta)
rN = size(db_psd,1);
n_seg = size(db_psd,2);
v = nan(rN,1);
snr = v;
noise_beam = v;
np = npts; % ex. 5 pts on each site of peak.
if plot_beamform_cur == 1
    figure
end
for ri = 1:rN  
%     y = movmean(db_psd(ri,:),fnpts);
    y = db_psd(ri,:);
    
    yp = y;
    yn = y;
    
    vp0 = (freq - f_bragg)*lambda_bragg; % velocities from positive Bragg peaks
    vn0 = (freq + f_bragg)*lambda_bragg; % velocities from negative Bragg peaks

    inoise1 = vp0 > vmax;
    noise1 = mean(y(inoise1),'omitnan');
    inoise2 = vn0 < -vmax;
    noise2 = mean(y(inoise2),'omitnan');
    noise_beam = (noise1 + noise2)/2;
    

    %%%%%%%%% for all frequencies within vmax, get average snr, do music db_min above this
    % pos
    fp = vp0 < vmax & vp0 > -vmax; % for positive Bragg peak in range +/- vmax
    yp(~fp) = min(yp)-1;
    [ypm,p] = max(yp);
    if isempty(p) || p == 1
        vp = nan;
    else
        fp = freq(p-np:p+np);
        yp1 = y(p-np:p+np);
        pp = polyfit(fp,yp1,2);
        fp1 = -pp(2)/(2*pp(1));
        if abs(freq(p) - fp1) > abs(fp(end)-fp(1)) % bad parabolic fit
            n   = 4; % 4th order weighting
            fp  = sum(fp.*yp1.^n)/sum(yp1.^n);  % peak frequency
        else
            fp = fp1; % polyfit
        end
        vp = (fp - f_bragg)*lambda_bragg;
    end
    
    % neg
    fn = vn0 < vmax & vn0 > -vmax; % for positive Bragg peak in range +/- vmax
    yn(~fn) = min(yn)-1;
    [ynm,n] = max(yn);
    if isempty(n) || n == 1
        vn = nan;
    else
        fn = freq(n-np:n+np);
        yn1 = y(n-np:n+np);
        pn = polyfit(fn,yn1,2);
        fn1 = -pn(2)/(2*pn(1)); % max
        if abs(freq(n) - fn1) > abs(fn(end)-fn(1)) % bad fit
            n   = 4; % 4th order weighting
            fn  = sum(fn.*yn1.^n)/sum(yn1.^n);  % peak frequency
        else
            fn = fn1; % use max bin
        end
        vn = (fn + f_bragg)*lambda_bragg;
    end
    
    % save
    if ypm >= ynm
        v(ri) = vp;
        f_max = fp;
        y_max = y(p);
        snr(ri) = y_max-noise1;
    else
        v(ri) = vn;
        f_max = fn;
        y_max = y(n);
        snr(ri) = y_max-noise1;
    end
    noise_beam(ri) = noise1;
    if plot_beamform_cur == 1
        hold off
        plot(freq,y)
        hold on

        y1 = y;
        y2 = y;
        y1(~inoise1) = nan;
        y2(~inoise2) = nan;
        plot(freq,y1)
        plot(freq,y2)
        plot(freq,noise_beam(ri)*ones(n_seg,1))
        plot(f_max,y_max,'ks')
        text(f_max,y_max,num2str(v(ri)))
        title(['range ' num2str(ri) ' angle' num2str(theta)])
        legend('signal','pos noise','neg noise','noise','max signal')
        drawnow

    end
end