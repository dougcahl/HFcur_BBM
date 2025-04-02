function radar_info_extra(f_radar)
lambda_radar = 299.8/f_radar;
lambda_bragg = lambda_radar/2;


% max peak velocity m/s +/- for the Bragg peak

% Bragg wavenumber, phase speed, and frequency
g = 9.81;
k_bragg = 2*pi/lambda_bragg;
v_bragg = sqrt(g/k_bragg);
w_bragg = sqrt(g*k_bragg);
f_bragg = w_bragg/2/pi;

v_bragg_harmonic = sqrt(2)*v_bragg;
v_bragg_harmonic_half = sqrt(1/2)*v_bragg;
v_bragg_harmonic3 = sqrt(3)*v_bragg;
v_bragg_harmonic_corner = (2^(3/4))*v_bragg;
v_bragg_harmonic_diff = v_bragg_harmonic-v_bragg;

f_bragg_harmonic = v_bragg_harmonic/lambda_bragg;
f_bragg_harmonic_half = v_bragg_harmonic_half/lambda_bragg;
f_bragg_harmonic3 = v_bragg_harmonic3/lambda_bragg;
f_bragg_harmonic_corner = v_bragg_harmonic_corner/lambda_bragg;

eff_depth = lambda_bragg/4/pi;
eff_depth_2x = 2*lambda_bragg/4/pi;
eff_depth_corner = sqrt(2)*lambda_bragg/4/pi;
eff_depth_3x = 3*lambda_bragg/4/pi;

disp(['radar frequency = ' num2str(f_radar) ' MHz'])
disp(' ')
disp(['lambda_bragg = ' num2str(lambda_bragg) ' m'])
disp(['lambda_bragg 2x harmonic = ' num2str(2*lambda_bragg) ' m'])
disp(['lambda_bragg 3x harmonic = ' num2str(3*lambda_bragg) ' m'])
disp(['lambda_bragg 1/2 harmonic = ' num2str(0.5*lambda_bragg) ' m'])
disp(['lambda_bragg corner = ' num2str(sqrt(2)*lambda_bragg) ' m'])
disp(' ')
disp(['bragg velocity = ' num2str(v_bragg) ' m/s'])
disp(['bragg 2x harmonic velocity = ' num2str(v_bragg_harmonic) ' m/s'])
disp(['difference = ' num2str(v_bragg_harmonic_diff) ' m/s'])
disp(['f_bragg = ' num2str(f_bragg)])
disp(['f 2x harmonic = ' num2str(f_bragg_harmonic)])
disp(['f 3x harmonic = ' num2str(f_bragg_harmonic3)])
disp(['f corner = ' num2str(f_bragg_harmonic_corner)])
disp(['f 1/2 harmonic = ' num2str(f_bragg_harmonic_half)])


disp(' ')
disp(['radar effective depth    = ' num2str(eff_depth) ' m'])
disp(['radar 2x effective depth = ' num2str(eff_depth_2x) ' m'])
disp(['radar 3x effective depth = ' num2str(eff_depth_3x) ' m'])
disp(['radar corner eff depth   = ' num2str(eff_depth_corner) ' m'])

% (sqrt(2)-1)*v_bragg = ~2.2 m/s at 8.3 MHz, so any peak this close to the Bragg peak
%could be the harmonic, i.e, not real data
