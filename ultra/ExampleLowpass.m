% EXAMPLE Design a lowpass filter with the following specs:
%
% Sampling frequency, rad/s.................: 20.0          
% Maximum passband ripple, Ap, in dB........: 3.0
% Minimum stoband attenuation, Aa, in dB....: 10.0         
% Passband edge wa rad/s....................: 4.0
% Stopband edge wp, rad/s...................: 6.0 

% Input parameters
clear all;
figure 
sampFreq = 20.0;
passRippleDB = 3.0;
stopRippleDB = 26;
passEdge = 4.0;
stopEdge = 6.0;

% Compute ripple for window design
stopRipple = 10^(-0.05*stopRippleDB);
passRipple = (10^(0.05*passRippleDB)-1) / (10^(0.05*passRippleDB)+1);
ripple = min(stopRipple,passRipple);

% The minimum stopband ripple allowable for the prediction equations is 25.4456 dB. 
% Anything less than this leads to a prediction for beta that is less than 1,
% however it is not possible to compute the ultraspherical window for beta < 1. 
% If someone asks for a larger ripple, he we set it to be 25.4456 dB.
minStopRippleDB = 25.4456; 
minStopRipple = 10^(-0.05*minStopRippleDB);
ripple = min(ripple,minStopRipple);

% Compute window design filter
[n,Wn,mu,beta,typ] = ultraord( [passEdge stopEdge], [1 0], [ripple ripple], sampFreq );
b = fir1(n, Wn, typ, ultra(n+1,mu,beta,'beta'), 'noscale');
freqz(b,1,512)

