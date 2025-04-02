function w = ultra(N, mu, par, partype)
%ULTRA Ultraspherical Window.
%   W = ULTRA(N,MU,BETA,'beta') returns the MU-valued N-point ultraspherical 
%   window with null-to-null width BETA times that of the rectangular window.
%   MU>-1.5, MU~=-1, and BETA>=1.
%
%   W = ULTRA(N,MU,XMU,'xmu') returns the MU- and XMU-valued N-point 
%   ultraspherical window. MU>-1.5, MU~=-1, and XMU>=1.
%
%   Notes:
%   The Dolph-Chebyshev and Saramaki windows are special cases of the 
%   ultraspherical window and can be obtained by letting MU = 0 and 1, 
%   respectively. Another special case is when MU = 0.5, which produces 
%   windows based on the Legendre polynomial.
%
%   Cautions:
%   Negative coefficients occasionally occur for extremely large/small values 
%   of BETA. Although this is not a problem, one can usually achieve a window 
%   with positive coefficients and similar frequency response by adjusting 
%   BETA appropriately.
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
%           Copyright 2005-2006 Stuart W.A. Bergen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input argument validation and some bound checking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error(nargchk(4,4,nargin));
narginchk(4,4);

% Check window length N
N = ceil(N);
if N < 1, error('The window length N must be >= 1.'); end
if N == 1, w=1; return; end % Trivial window

% Check parameter MU
if mu <= -1.5 | mu == -1, error('MU must be > -1.5 and ~= -1.'); end

% Check parameter BETA/XMU
partypelength = length(partype);
if partypelength == 4,
    if partype == 'beta',
        if par < 1, error('BETA must be >= 1.'); end
        xmu = largestZero(N-1,mu)/cos(pi*par/N); % calculate xmu from beta
    else
        error('PARTYPE must be ''beta'' or ''xmu''.');
    end
elseif partypelength == 3,
    if partype == 'xmu',
        if par < 1, error('XMU must be >= 1.'); end
        xmu = par; % par supplied is xmu
    else
        error('PARTYPE must be ''beta'' or ''xmu''.');
    end
else
    error('PARTYPE must be ''beta'' or ''xmu''.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate window coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate constants.
isodd=mod(N,2); if isodd, A=(N-1)/2; else A=N/2-1; end
p=N-1; Cval=1-xmu^-2;
if mu==0, B=xmu^p; else B=mu*xmu^p; end % Limiting value for Chebyshev polynomial.

% Working vector initalizations.
v0(1,A+1)=0; v1(1,A+1)=0; v2(1,A+1)=0; C(1,A+1)=0; 
ww(1,A+1)=0; w(1:N)=0;

% Generate fixed v0(n+1) and C(n+1)^n for n=0,1,...,A. 
a0=mu+p-1; b0=p-1; v0(1)=1; C(1)=1; 
for g=0:b0-1, v0(1)=(a0-g)/(b0-g)*v0(1); end % v0(1)=binomial(a0,b0).
for n=1:A, nl=n+1; v0(nl)=(b0-n+1)/(a0-n+1)*v0(nl-1); C(nl)=Cval*C(nl-1); end

% Calculate 'half' the window coefficients ww(n+1) for n=0,1,...,A.
for n=0:A, nl=n+1;
	% Generate v1(g+1) for g=n,n-1,...,0 and v2(g+1) for g=0,1,...,n.
	v1(nl)=1; v2(1)=1; alpha2=p-n;
	for g=1:n, gl=g+1; m=n-g; ml=m+1;
        v1(ml)=(mu+m)/(n-m)*v1(ml+1); v2(gl)=(alpha2-g+1)/g*v2(gl-1);
	end
	sum1=sum((v1(1:nl).*v2(1:nl)).*C(1:nl)); % Perform sumation.
	ww(nl)=B/alpha2*v0(nl)*sum1; % Window coefficients.
end

% Obtain symmetrical normalized window coefficients.
ww=ww/ww(end);
if isodd, w=[ww ww(end-1:-1:1)]';
else w=[ww fliplr(ww)]'; end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LARGESTZERO computes the largest zero of the ultraspherical polynomial 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xstar = largestZero(n,alpha)

% Check bounds for parameters.
if (floor(n)~=n) | (n<=1), error('N should be a positive integer > 1'); xstar=NAN; return; end
if alpha<=-1.5, error('ALPHA should be >= -1.5'); xstar=NAN; return;
elseif alpha==-1 error('ALPHA cannot be -1'); xstar=NAN; return; end

% Check for special instances of the ultraspherical polynomial where 
% analytical expressions for the zeros exist.
if alpha==0, xstar=cos(pi/(2*n)); return; % Chebyshev polynomial of the first kind, Tn(x).
elseif alpha==1, xstar=cos(pi/(n+1)); return; % Chebyshev polynomial of the second kind, Un(x).
elseif alpha==-0.5, xstar=1; return; end % Always 1.  

% Step 1
epsilon=10^-6; ktol=20; C(1,n)=0; y(1,ktol+1)=0;
y(1)=sqrt(n*n+2*n*alpha-2*alpha-1)/(n+alpha); % Upper bound

% Step 2
for k=1:ktol
	x=y(k);	C(1)=2*alpha.*x; C(2)=-alpha+2*alpha*(1+alpha).*(x.^2);
	for ii=3:n,
	    C(ii)=2*(ii+alpha-1)/(ii).*x.*C(ii-1)-(ii+2*alpha-2)/(ii).*C(ii-2); 
    end 
	den=1/(1-x^2)*((2*alpha+n-1)*C(n-1)-(n*x)*C(n)); % derivative.
	y(k+1)=y(k)-C(n)/den; % Newton-Raphson iteration.
	
% Step 3
	if abs(y(k+1)-y(k))<epsilon, xstar=y(k+1); break; end
end
if k==ktol, error('Algorithm did not converge within ten iterations - stopping.'); end
return
%--------------------------EOF-------------------------------