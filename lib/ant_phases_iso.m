function [zphi,theta1,theta_radar] = ant_phases_iso(beam_lim,fradar,lon,lat,dtheta,r,s_plot,s_plotname,theta_guess)

lambda_radar = 299.8/fradar; % meters

theta1 = beam_lim(1):dtheta:beam_lim(2);


% antenna positions
mean_lon = mean(lon);
mean_lat = mean(lat);
[x_radar,y_radar] = geog2utm(lon,lat,mean_lon,mean_lat); % zero at midpoint of antennas
x_radar = x_radar*1000;
y_radar = y_radar*1000;

% get angles of beam directions
p_radar = polyfit(x_radar,y_radar,1);
slope_radar = p_radar(1);
theta_guess_math = 90 - theta_guess; % -270 : 90

theta_radar = atan2d(slope_radar,1) - 90; % -270:90
dth = abs(theta_radar - theta_guess_math);
if dth > 90
    theta_radar = atan2d(slope_radar,1) + 90; % -90:270
    dth = abs(theta_radar - theta_guess_math);
    if dth > 180
        dth = abs(dth - 360);
    end
    if dth > 90
        error('bad guess of angle')
    end
end
theta = theta1 + theta_radar;

if s_plot == 1
    figure;
    plot(x_radar,y_radar,'-k o')
    hold on
    plot(x_radar(1),y_radar(1),'r*')
    plot(0,0,'kx')
    axis equal
    quiver(0,0,50*cosd(theta_radar),50*sind(theta_radar),0)
    text(50*cosd(theta_radar),50*sind(theta_radar),['\theta = ' num2str(theta_radar) '^\circ'])
    xlabel('m')
    ylabel('m')
    box off
    grid on
    legend('antennas','antenna 1','center of array')
    title(s_plotname)
    text(1.04,1,'N','units','normalized')
    annotation('arrow',[.94 .94],[.8 .9])
    drawnow
end

%%%%%%%%%%%%%%%%%%%%


zphi = nan(length(x_radar),length(theta),length(r));
for j = 1:length(r)
    dist = r(j); % in meters
    x = dist*cosd(theta);
    y = dist*sind(theta);
    
    %%
    zphi1 = nan(length(x_radar),length(theta));
    for i=1:length(theta)
        d = sqrt( (x(i)-x_radar).^2 + (y(i)-y_radar).^2);
        dd = d - d(1);
        
        dphi = -2*pi*dd/lambda_radar;
        zphi1(:,i) = exp(1i.*dphi);
    end
    zphi(:,:,j) = zphi1;
end

end