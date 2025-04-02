function radar_info(f_radar)
lambda_radar = 299.8/f_radar;
lambda_bragg = lambda_radar/2;


% max peak velocity m/s +/- for the Bragg peak

% Bragg wavenumber, phase speed, and frequency
g = 9.81;
k_bragg = 2*pi/lambda_bragg;
v_bragg = sqrt(g/k_bragg);
w_bragg = sqrt(g*k_bragg);
f_bragg = w_bragg/2/pi;
eff_depth = lambda_bragg/4/pi;

disp(['Radar frequency  = ' num2str(f_radar) ' MHz'])
disp(['Radar wavelength = ' num2str(lambda_radar) ' m'])
disp(['Bragg wavelength = ' num2str(lambda_bragg) ' m'])
disp(['Bragg velocity   = ' num2str(v_bragg) ' m/s'])
disp(['Bragg frequency  = ' num2str(f_bragg) ' Hz'])
disp(['Effective depth  = ' num2str(eff_depth) ' m'])
