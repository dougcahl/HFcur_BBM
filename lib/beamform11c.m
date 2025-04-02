function zbeam= beamform11c(zz,as,w_beam)
% function zbeam= beamform11c(zz,as,w_beam)
%
% computes beamform from I,Q data from WERA .SORT file using
% complex antenna pattern zph and a beamforming window w_beam
%
% input
% zz - covariance matrices of FFT of RX signal, .SORT from WERA
% a(antN,thN,rN) - complex antenna pattern
% w_beam(antN) - beamforming window
%
%
% output
% zbeam(rN,chirpN,thN) - beamformed amplitude
%
% DLC 2023
%

rN = size(zz,1);    % range cells
fN = size(zz,2);    % freq bins 
thN = size(as,2);   % th bins

% u_beam = sum(w_beam);

% preallocate
zbeam = nan(rN,fN,thN);

for ic = 1:fN
    for ir = 1:rN
        R = zz{ir,ic};
        a = as(:,:,ir).*w_beam;
        zbeam(ir,ic,:) = abs(diag(a'*R*a./(a'*a)));
%         for ith = 1:thN
%             a = as(:,ith,ir).*w_beam;
%             zbeam(ir,ic,ith) = abs(a'*R*a/(a'*a));
%         end
    end
end