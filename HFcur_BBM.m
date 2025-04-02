%% HFcur_BBM.m
% A MATLAB package for calculating surface currents from HF
% radars using Beamforming, Beamscan (direction finding using beamforming) 
% and MUSIC for a linear array (v0.0aa). 
% Zenodo. https://doi.org/10.5281/zenodo.7231459
%
% Method described in Cahl, D., G. Voulgaris, and L. Leonard, 2023:
% A Comparison of Beamforming and Direction Finding Algorithms (Beamscan 
% and MUSIC) on a Linear Array HF Radar in a Medium to Low Wave Energy 
% Environment. J. Atmos. Oceanic Technol., 40, 191–218,
% https://doi.org/10.1175/JTECH-D-22-0005.1. 
%
% Douglas Cahl and George Voulgaris 2023
% University of South Carolina
%
%% Toolboxes used
% matWERA 
% George Voulgaris, 2019. 
% matWERA: A MATLAB package for reading binary data files recorded by a 
% Wellen-Radar (WERA) HF Radar. 
% Zenodo. http://doi.org/10.5281/zenodo.3570891.
%
% Ultraspherical window 
%  Stuart W.A. Bergen, 12-8-2005. sbergen@ece.uvic.ca
%  This file is an adapted version of KAISERORD written by 
%  J.H. McClellan, 10-28-1991, Revision: 1.10 Date: 2002/03/28, Mathworks.
%
%
% please set parameters in this file to run Beamforming, Beamscan and MUSIC
% 

%% add toolbox paths
addpath('matWERA') % for reading WERA .SORT files
addpath('ultra')   % for ultraspherical FFT window
addpath('lib')     % for subroutines
% addpath('MATLAB/codegen/mex/hf_beam_max_current_upd1')
% addpath('MATLAB/codegen/mex/beamform1')


%% setup
theta_guess = 180; % guess of boresight direction degrees North (0-360)
% calculates exact slope from antenna locations file

beam_lim = [-60 60]; % degrees from boresight to analyse data

% beamform
run_beamforming         = 1; % =1 run beamforming, =0 don't run
save_data_beamform      = 1; % save beamformed spectra (creates large file)
save_data_beamform_cur  = 1; % save Beamformed currents 

% beamscan
run_beamscan    = 1; % run beamscan
save_beamscan   = 1; % save Beamscan currents

% music
run_music       = 1; % run MUSIC
save_music      = 1; % save MUSIC currents

%% overall plots
ant_plot            = 1; % plot antenna locations for each SORT file and 
% boresight direction in mathematical degrees
fnpts               = 9; % number of points for movmean smoothed FFT plots 

plot_DOA_analysis   = 2; % plot comparing results from Beamforming,
% Beamscan and MUSIC for each range. =1 plots of DOAs, =2 plots of DOAs
% with signal strength. I recommend using 2 for analysis.

%% radar information 
site    = 'gtn';                    % site name
% site    = 'csw';                    % site name
sortdir   = ['data/' site];         % sort files directory
savedir = ['data/processed/' site]; % save directory

ddsclk = 90.742153846154*2*10^6; % find in WERA.mes (needed for range processing)
% This is likely common to systems with the 90 MHz clock?

% antenna locations file (from WERA PC) for ant_lons, ant_lats
% ant_pos_file = 'data/antpos_CSW.asc';
ant_pos_file    = 'data/antpos_GTNv5.asc';
ant             = importdata(ant_pos_file,' ',1); % it has 1 header line
ant_lons        = ant.data(:,3); 
ant_lats        = ant.data(:,2);
antN            = length(ant_lons); % number of antennas ex. 8,12, etc

vmax = 2.0; % max velocity m/s +/- that is possible in your radar coverage area 
% in the future vmax this should be adjustable for different locations

% f_radar does not have to be exact as it will read it from the SORT file
f_radar = 8.3; % radar operating frequency in MHz 
radar_info(f_radar); % display radar info
% radar_info_extra(f_radar); % display more radar info

%% load measured antenna pattern (not implelented)
% load(['ant_pattern_' site '.mat']); % 
% measured_pattern = 1;

%% FFT window
window_FFT = 1; % 1 = Blackman Harris

%% Beam forming setup 
window_BF       = 1;    % 1 for Hamming, 2 for ultraspherical
n_BF            = 512;  % FFT sample length
p_overlap_BF    = 50;   % FFT percent overlap 

% beamforming options
npts            = 5; % number of points on each side of peak for frequency peak calc
dtheta_beam     = 1;    % save/process data every dtheta_beam degrees

% plots
plot_doppler_Bform      = 1; % plot of antenna and avg spectra for Beamforming
plot_doppler_Bform_ants = 1; % plots range Doppler for each antenna for Beamformign
plot_beamform_cur       = 0; % plot for each current fit (many plots!)

%% Beamscan setup
window_BS       = 1;    % 1 for Hamming, 2 for ultraspherical

% covariance matrix setup 
% this is the highest velocity resolution Beamscan for this data set
n_BS            = 2048;
p_overlap_BS    = 0;

% example alternate 1, lower resolution but maybe better signal quality
% n_BS = 1024;
% p_overlap_BS = 50;

% example alternate 2, even lower resolution but maybe even better signal
% quality
% n_BS = 512;
% p_overlap_BS = 50;

db_min_BS      = 5; % min dB for Bragg peak
bscan_usesmoothed_spec = 0; % use a fnpts movemean for SNR identification of freq bins for Beamscan

% plots
plot_doppler_beamscan           = 1; % plot of antenna and avg spectra for Beamscan
plot_db_peaks_BS                = 0; % plot the peaks and corresponding frequency bins that will be used to calculate surface currents for each range
plot_BS_ind_spec                = 0; % plot the Beamscan spectrum for each range/freqnecy to calculate a DOA
plot_beamscan_DOAs_per_range    = 0; % plot of all DOAs for each range

%% MUSIC setup
% MUSIC, modified MUSIC (MMUSIC) or improved modified MUSIC (IMMUSIC)
MUSIC_method = 1; % 1 = MUSIC,  2 = MMUSIC,  3 = IMMUSIC 
% 2 and 3 are only partially implemented, the matrices are calculated but
% the full implementation into this script is not, see MUSIC.m for the MMUSIC
% and IMMUSIC sections that are commented out

% (important tuning parameter)
nmax        = 4; % max number of DOAs, must be less than antenna number, antN!

% Kirincich peakfinder threshold (important tuning parameter)
peak_thresh = 2.0; % 2.0; % increase for more antennas, closer to 1/2 wavelength spacing. T
% decrease for less antennas, closely spaced antennas, nonlinear arrays
% this is a critical tuning parameter for MUSIC to work for HF radar data!
% Spend time adjusting this

db_min_MUSIC      = 5; % min dB for Bragg peak
music_usesmoothed_spec = 0; % use a fnpts movemean for SNR identification of freq bins for MUSIC


% covariance matrix setup
n_MU            = 1024; % segment length of FFT 
p_overlap_MU    = 92; % percent overlap

% plots 
plot_doppler_MUSIC          = 1; % plot of antenna and avg spectra for MUSIC 
plot_db_peaks_MUSIC         = 0; % plot the peaks and corresponding frequency bins that will be used to calculate surface currents for each range
plot_MUSIC_DOAs             = 0; % plot the MUSIC psuedospectrum for each range/freqnecy to calculate a DOAs
plot_MUSIC_DOAs_per_range   = 0; % plot of all DOAs for each range


%% process Doppler spectra
HFcur_BBM_sub;

%% processing currents to a filled in current map 
% I was using this code for single locations and the full current maps are
% a work in progress currently.

% after this you want to interpolate onto all the radar current
% measurement locations. 

% For beamforming, there is a results at almost all
% locations but maybe a minimum SNR will be helpful to reduce outliers
% and/or some weighted interpolation function over space to remove outliers

% For Beamscan and MUSIC, either a signal (for MUSIC) or SNR (for beamscan)
% weighted interpolation function over space or a signal/SNR threshold to
% remove outliers will be necessary. Beamscan and MUSIC produce far fewer
% DOAs than beamforming so interpolation in space (common) and/or time
% (uncommon) is required for a full current map. 

% There are several other methods for outlier removal in the literature

