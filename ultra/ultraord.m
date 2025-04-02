function [N, Wn, mu, beta, typ] = ultraord(fcuts, mags, devs, fsamp, cellflag)
%ULTRAORD FIR order estimator (lowpass, highpass, bandpass, multiband).
%   [N,Wn,MU,BETA,TYPE] = ULTRAORD(F,A,DEV,Fs) is the approximate order N, 
%   normalized frequency band edges Wn, ultraspherical window parameters MU
%   and BETA, and filter type TYPE to be used by the FIR1 function:
%      B = FIR1(N, Wn, TYPE, ultra( N+1,MU,BETA,'beta' ), 'noscale' )
%
%   The resulting filter will approximately meet the specifications given
%   by the input parameters F, A, and DEV.
%
%   F is a vector of band edge frequencies in Hz, in ascending order between 0
%   and half the sampling frequency Fs.  A is a vector of 0s and 1s specifying 
%   the desired function's amplitude on the bands defined by F. The length of F
%   is twice the length of A, minus 2 (it must therefore be even).  The first 
%   frequency band is assumed to start at zero, and the last one always ends 
%   at Fs/2.
%
%   DEV is a vector of maximum deviations or ripples allowable for each band.
%   The smallest deviation specified (MIN(DEV)) is used for both the passband
%   and the stopband.
%
%   Fs is the sampling frequency (which defaults to 2 if you leave it off).
%
%   C = ULTRAORD(F,A,DEV,Fs,'cell') is a cell-array whose elements are the
%   parameters to FIR1.
%
%   EXAMPLE
%      Design a lowpass filter with a passband edge of 1500Hz, a 
%      stopband edge of 2000Hz, passband ripple of 0.01, stopband ripple 
%      of 0.1, and a sampling frequency of 8000Hz:
%
%      [n,Wn,mu,beta,typ] = ultraord( [1500 2000], [1 0], [0.01 0.1], 8000 );
%      b = fir1(n, Wn, typ, ultra(n+1,mu,beta,'beta'), 'noscale');
%   
%      This is equivalent to
%      c = ultraord( [1500 2000], [1 0], [0.01 0.1], 8000, 'cell' );
%      b = fir1(c{:});
%
%   CAUTION 1: The order N is just an estimate. If the filter does not
%   meet the original specifications, a higher order such as N+1, N+2, etc. 
%   will; if the filter exceeds the specs, a slightly lower order one may work.
%   CAUTION 2: Results are inaccurate if cutoff frequencies are near zero
%   frequency or the Nyquist frequency; or if the devs are large (10%).
%
%   References:
%   [1] S.W.A. Bergen, "Design of the Ultraspherical Window Function and Its 
%       Applications," PhD Dissertation, University of Victoria, Sept. 2005.
%   [2] S.W.A. Bergen and A. Antoniou, "Design of Ultraspherical Window 
%       Functions with Prescribed Spectral Characteristics," EURASIP Journal 
%       on Applied Signal Processing, vol. 2004, no. 13, pp. 2053-2065, 2004.
%   [3] S.W.A. Bergen and A. Antoniou, "Design of Nonrecursive Digital Filters 
%       Using the Ultraspherical Window Function," EURASIP Journal on Applied 
%       Signal Processing, vol. 2005, no. 12, pp. 1910-1922, 2005.

%   Author: Stuart W.A. Bergen, 12-8-2005. sbergen@ece.uvic.ca
%           This file is an adapted version of KAISERORD written by 
%           J.H. McClellan, 10-28-1991, Revision: 1.10 Date: 2002/03/28, 
%           Mathworks.

error(nargchk(3,5,nargin))
if nargin < 4 | isempty(fsamp),
    fsamp = 2;
end
if max(fcuts) >= fsamp/2,
    error('Band edge frequencies must between between 0 and half the sampling frequency Fs.');
end

fcuts = fcuts/fsamp;       %  NORMALIZE to sampling frequency

% Turn vectors into column vectors
fcuts = fcuts(:);
mags = mags(:);
devs = devs(:);

mf = size(fcuts,1);
nbands = size(mags,1);

if size(mags,1) ~= size(devs,1)
    error('Requires A and DEV to be the same length.');
end
if( min(abs(mags)) )
   error('Stopbands must be zero.');
end
dmags = abs(diff(mags));
if( any(dmags~=dmags(1)) )
    error('All passbands must have same height.');
end
if( any(diff(fcuts)<0) )
    error('Bandedges must be strictly increasing.');
end


if mf ~= 2*(nbands-1)
    error('Length of F must be 2*length(A)-2.');
end

zz = mags==0;             % find stopbands
devs = devs./(zz+mags);   % divide delta by mag to get relative deviation

% Determine the smallest width transition zone
% Separate the passband and stopband edges
%
f1 = fcuts(1:2:(mf-1));
f2 = fcuts(2:2:mf);
[df,n] = min(f2-f1);

%=== LOWPASS case: Use formula (ref: Herrmann, Rabiner, Chan)
%
if( nbands==2 )
     [L,mu,beta] = ultralpord( f1(n), f2(n), devs(1), devs(2));

%=== BANDPASS case:
%    - try different lowpasses and take the WORST one that
%        goes thru the BP specs; try both transition widths
%    - will also do the bandreject case
%    - does the multi-band case, one bandpass at a time.
%    
else
  L = 0;  mu = 0; beta = 0;
  for i=2:nbands-1,
    [L1,mu1,beta1] = ultralpord( f1(i-1), f2(i-1), devs(i),   devs(i-1) );
    [L2,mu2,beta2] = ultralpord( f1(i),   f2(i),   devs(i),   devs(i+1) );
    if( L1>L )
        mu = mu1; beta = beta1;  L = L1;   end
    if( L2>L )
        mu = mu2; beta = beta2;  L = L2;   end
  end
end

N = ceil( L ) - 1;   % need order, not length, for Filter design

%=== Make the MATLAB compatible specs for FIR1
%
Wn = 2*(f1+f2)/2;    %-- use mid-frequency; multiply by 2 for MATLAB
typ = 'low';
if( nbands==2 & mags(1)==0 )
  typ='high';
elseif( nbands==3 & mags(2)==0 )
  typ='stop';
elseif( nbands>=3 & mags(1)==0 )  
  typ='DC-0';                    
elseif( nbands>=3 & mags(1)==1 ) 
  typ='DC-1';                   
end

% If order is odd, and gain is not zero at nyquist, increase the order by one.
if rem(N,2) & mags(end)~=0,
    N = N + 1;
end

if nargout == 1 & nargin == 5
  N = {N, Wn, typ, ultra(N+1,mu,beta,'beta'), 'noscale'};
end
return
%%%% ---- end of ultraord


function [L,mu,beta] = ultralpord(freq1, freq2, delta1, delta2 )
%ULTRALPORD FIR lowpass filter Length estimator
%
%   [L,mu,beta] = ultralpord(freq1, freq2, dev1, dev2)
%
%   input:
%     freq1: passband cutoff freq (NORMALIZED)
%     freq2: stopband cutoff freq (NORMALIZED)
%      dev1: passband ripple (DESIRED)
%      dev2: stopband attenuation (not in dB)
%
%   outputs:
%      L = filter Length (# of samples)   **NOT the order N, which is N = L-1
%   beta =  parameter for the ultrasperical window
%
%   NOTE: Will also work for highpass filters (i.e., f1 > f2)
% 	      Will not work well if transition zone is near f = 0, or
%         near f = fs/2

% 	Author(s): Stuart W.A. Bergen, 12-8-2005.
%              This file is an adapted version of KAISLPORD written by
%              J.H. McClellan, 8-28-95
	
%   References:
%     [1] Rabiner & Gold, Theory and Applications of DSP, pp. 156-7.     

delta = min( [delta1,delta2] );
atten = -20*log10( delta );
if atten<=80, D = 4.6447e-05*atten.^2 + 0.062161*atten -0.48184;
else,  D = 1.7098e-06*atten.^2 + 0.070894*atten -0.89368; end
df = abs(freq2 - freq1);
L = D/df + 1;
mu =-1.7211e-05*atten.^2 + 6.7213e-03*atten + 0.18974;
if atten<60, beta = 0.00004024*atten^2 + 0.02423*atten + 0.3574;
elseif atten>=60 & atten<120, beta = 0.00007303*atten^2 + 0.02079*atten + 0.4447; 
elseif atten>=120, beta = 6.7326e-06*atten^2 + 0.033374*atten -0.1192; end
return
%--------------------------EOF-------------------------------
