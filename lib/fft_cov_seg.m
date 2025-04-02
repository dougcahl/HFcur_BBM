function [freq,db_psd,zz_mat1,zz_mat2,zz_mat3] = fft_cov_seg(I,Q,w,fs,n_seg,p_overlap)
%#codegen
%
antN = size(I,1);
rN = size(I,2);
% chirpN = size(I,3);
u = 1/length(w)*sum(w.^2);
freq = 0;
% FFT calcs
db_psd = nan(antN,rN,n_seg);
% zz = complex(db_psd);
zz_mat = cell(rN,n_seg);
zz_mat1 = coder.nullcopy(zz_mat);
zz_mat2 = coder.nullcopy(zz_mat);
zz_mat3 = coder.nullcopy(zz_mat);
data = squeeze(I(1,1,:)) + 1i*squeeze(Q(1,1,:));
y   = make_segments_overlap(data,n_seg,p_overlap);%,0);
m = size(y,2); % number of segments
% zz_mat1_segs = complex(nan(rN,antN,ys2,n_seg + 1));
for ir = 1:rN
    zz1_00 = complex(nan(antN,m,n_seg));
    % phase_00 = nan(antN,m,n_seg);
    % energy_00 = nan(antN,m,n_seg);
    for in = 1:antN
        data = squeeze(I(in,ir,:)) + 1i*squeeze(Q(in,ir,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%% FFT averaging in 512 or 1024 segments
        df  = fs/n_seg;             % Spectral resolution
        %             bandwidth = (7/6)*df;   % effective bandwidth
        y   = make_segments_overlap(data,n_seg,p_overlap);%,0);
        %             n_samples = size(y,2);
        m   = size(y,2);    % number of segments
        % energy  = nan(m,n_seg);
        zz1 = complex(nan(m,n_seg));
        % phase = nan(size(y,2),n_seg);
        % amp = nan(size(y,2),n_seg);
        for ii = 1:m
            y1 = y(:,ii).*w;         % multiply the segment by window
            % [freq,energy(ii,:),amp(ii,:),~,phase(ii,:),zz1(ii,:)] = fft_doug_phase_2(y1, fs);
            [zz1(ii,:),freq,df] = fft_doug_phase_2(y1, fs); % FFT
        end
        zz1_00(in,:,:) = zz1;
        db_psd(in,ir,:) = 10*log10(mean(abs(zz1.^2)/u/df));
     
    end
    % now we can make matrices using all 12 antennas for avg FFTs
    for izz = 1:n_seg
        % reference antenna 1 each time for cov matrix?
        zzt1 = zz1_00(:,:,izz);
        zzt2 = zzt1*zzt1';
        zzt2 = sqrt(1/u)*zzt2/m;
        R = zzt2;
        J=fliplr(eye(antN));
        Ry = J*conj(R)*J;
%         R=(R+J*conj(R)*J)/2; % other cov matrix
        zz_mat1{ir,izz} = R; % improved MUSIC cov matrix
        zz_mat2{ir,izz} = (R+Ry)/2; % improved MMUSIC cov matrix
        zz_mat3{ir,izz} = [R Ry]; % improved IMMUSIC cov matrix
    end
end

