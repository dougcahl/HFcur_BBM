function [r,f_offset] = wera_offsets(R,WERA,ddsclk)
% WERA frequency and range offset

dr = mode(diff(R)); % range bin size
foff = -1*WERA.RXOFFSET*ddsclk/(2^48);
if WERA.T_chirp > 10^4
    rf = 1/(WERA.T_chirp/10^6);
else
    rf = 1/(WERA.T_chirp);
end

fs = 1/WERA.RATE;
roff = floor(foff/rf);
doff = foff-roff*rf;

if doff > rf/2
    roff = roff + 1;
    doff = foff-roff*rf;
end

% frequency offset
if abs(abs(doff) - fs) < .01
    f_offset = 0;
else
    f_offset = -doff; % Doppler freq = freq - f_offset;
end

% range
R = R + roff*dr;
% convert meters
r = R*1000;

