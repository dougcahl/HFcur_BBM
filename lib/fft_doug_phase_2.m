function [Y,freq,df] = fft_doug_phase_2(data, Fs)
% positive and negative frequencies
% function [Y,freq,df] = fft_doug_phase_2(data, Fs)
%
%   Inputs
%           data    = your data
%           Fs      = sampling frequency
%
%   Outputs
%           freq    = frequencies
%           df      = Spectral resolution
%           Y       = complex amplitudes fftshifted to match frequencies
%
% DLC 2014
%

L           = length(data);             % Length of signal
df          = Fs/L;                     % Spectral resolution
freq        = (-L/2:L/2-1)*df;          % Frequencies

Y           = fftshift(fft(data))/L;    % Calculate fft Fourier Coefficients

end
