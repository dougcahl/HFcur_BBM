function [theta00,signal1] = beamscan_single_bragg_peakc(ant_pattern_phi,freq00,zz00,theta1,plot_BS_ind_spec,tname,w_beam)
%#codegen
theta00 = nan(length(freq00),1);
signal1 = nan(length(freq00),1);
if plot_BS_ind_spec == 1
    figure
end
for i_00 = 1:length(freq00)
    % Beamscan: scan the beam over theta for each frequency
    zz = cell2mat(zz00(:,i_00)); % complex value in each antenna
%     I1 = real(zz);
%     Q1 = imag(zz);
%     pt = nan(size(theta1));
%     for ien = 1:length(theta1)
%         zp = ant_pattern_phi(:,ien); % a(theta)
%         Ex = sum(I1.*real(zp).*w_beam+Q1.*imag(zp).*w_beam);
%         Ey = sum(-I1.*imag(zp).*w_beam+Q1.*real(zp).*w_beam);
%         Ecom = squeeze(Ex + 1i*Ey);
%         pt(ien) = abs(Ecom).^2; 
%     end
    a1 = w_beam.*zz;
    pt = abs(a1'*ant_pattern_phi).^2;
    
    w_sum = (sum(w_beam).^2);
    ant_sum = sum(abs(sum(ant_pattern_phi,1).^2));
    pt = pt/w_sum/ant_sum;
    pt = sum(pt,1);
    
    % get peak max
    [ptm,pti] = max(pt);
    signal1(i_00) = ptm; % snr
    theta00(i_00) = theta1(pti);

    if plot_BS_ind_spec == 1
        hold off
        plot(theta1,pt)
        hold on
        plot(theta00(i_00),signal1(i_00),'rs')
        title(tname)
        drawnow
    end
end % for each frequency peak

