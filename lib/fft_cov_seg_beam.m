function [freq,zz] = fft_cov_seg_beam(I,Q,w,fs,n_seg,p_overlap)
% fft with ant1 reference and average of segments out
%#codegen
% [db,zz,zzmat] = function fft_cov(data,wn3,n_seg,p_overlap)
%
antN = size(I,1);
rN = size(I,2);
chirpN = size(I,3);
u = 1/length(w)*sum(w.^2);
freq = 0;
% FFT calcs
% db_psd = nan(antN,rN,n_seg);
% zz = complex(db_psd);
zz = zeros(antN,size(I,2),n_seg);
data = squeeze(I(1,1,:)) + 1i*squeeze(Q(1,1,:));
y   = make_segments_overlap(data,n_seg,p_overlap,0);
ys2 = size(y,2); % number of segments
% zz_mat1_segs = complex(nan(rN,antN,ys2,n_seg + 1));
for ir = 1:rN
    zz1_00 = complex(nan(antN,ys2,n_seg));
    phase_00 = nan(antN,ys2,n_seg);
    energy_00 = nan(antN,ys2,n_seg);
    for in = 1:antN
        data = squeeze(I(in,ir,:)) + 1i*squeeze(Q(in,ir,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%% FFT averaging in 512 or 1024 segments
        df  = fs/n_seg;             % Spectral resolution
        %             bandwidth = (7/6)*df;   % effective bandwidth
        y   = make_segments_overlap(data,n_seg,p_overlap,0);
        %             n_samples = size(y,2);
        m   = size(y,2);    % number of segments
        energy  = nan(m,n_seg);
        zz1 = complex(nan(m,n_seg));
        phase = nan(size(y,2),n_seg);
        amp = nan(size(y,2),n_seg);
        for ii = 1:m
            y1 = y(:,ii).*w;         % multiply the segment by window
            [freq,energy(ii,:),amp(ii,:),~,phase(ii,:),zz1(ii,:)] = fft_doug_phase_2(y1, fs);
        end
        zz1_00(in,:,:) = zz1;
        phase_00(in,:,:) = phase;
        energy_00(in,:,:) = energy;
        %             L = length(data);
        %             EDF         = 2.8*L*(7/6)/(N);    % equivalent degrees of freedom
        %             if m > 1
        energy_1d      = (1/u)*mean(energy);
        zz1_1d      = sqrt(1/u)*mean(zz1,1); %%%%%%%%%%%%%%%%%%%%%%
        zz1_1d = zz1_1d(:);
        %             else
        %                 energy_1d      = (1/u)*energy;
        %                 zz1_1d      = sqrt(1/u)*zz1;
        %             end
        psd         = energy_1d./df;
        %             amp         = sqrt(energy);
%         db_psd1 = 10*log10(psd);
%         db_psd(in,ir,:) = db_psd1(1:end);
%         zz(in,ir,:) = zz1_1d(1:end);
    end
    % now we can make matrices using all 12 antennas for avg FFTs
    for izz = 1:m
        aa1 = squeeze(zz1_00(:,izz,:)); % segment FFT vals
        a1 = aa1(1,:); % reference antenna 1
        p1 = atan2(imag(a1),real(a1)); % antenna 1 phase
        zz1_00(:,izz,:) = aa1.*exp(-1i*p1); % phase shift to ant1 = 0 phase
    end
    zz(:,ir,:) = squeeze(mean(zz1_00,2))*sqrt(1/u)/m;

end





















