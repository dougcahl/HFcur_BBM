function [freq,zz,db_psd,db_psd_ant,df] = fft_beam(I,Q,fs,w,n_seg,p_overlap)
%%%%%%%%%%%%%%%%%%% need to make this just spit out ffts, then average over
%%%%%%%%%%%%%%%%%%% the covariance matrices for averaging instead

% Blackman Harris window for segment FFTs
% wn = window(@blackmanharris,n_seg); 

antN = size(I,1);       % number of antennas
rN = size(I,2);         % range cells
chirpN = size(I,3);     % chirps

% FFT calcs
u = 1/length(w)*sum(w.^2); % normalizing factor

% preallocate
db_psd = nan(rN,n_seg);
zz = cell(rN,n_seg);
db_psd_ant = nan(antN,rN,n_seg);
% zz2 = cell(rN,n_seg);
% zz3 = cell(rN,n_seg);

% calculate number of segments
dx  = round(n_seg*(1-p_overlap/100));  
m   = floor((chirpN-n_seg)/dx) + 1;     % number of segments

for ir = 1:rN
    zz1_00 = complex(nan(antN,m,n_seg));
    for in = 1:antN
        zz1 = complex(nan(m,n_seg));
        data = squeeze(I(in,ir,:)) + 1i*squeeze(Q(in,ir,:));
        y   = make_segments_overlap(data,n_seg,p_overlap); 
        for ii = 1:m
            y1 = y(:,ii).*w;     % multiply the segment by window
            [zz1(ii,:),freq,df] = fft_doug_phase_2(y1, fs); % FFT
        end
        zz1_00(in,:,:) = zz1; % save FFT
        if m > 1
            db_psd_ant(in,ir,:) = 10*log10(mean(abs(zz1.^2)/u/df));
        else
            db_psd_ant(in,ir,:) = 10*log10(abs(zz1.^2)/u/df);
        end
    end

    % now we can make matrices using all 12 antennas for avg FFTs
    % zz1_00_mat = complex(nan(antN,antN,m));
    for izz = 1:n_seg
        zzt1 = zz1_00(:,:,izz);
        R = zzt1*zzt1'/u/m;      
%         J=fliplr(eye(antN));
%         Ry = J*conj(R)*J;
        zz{ir,izz} = R; % cov matrix for Beamscan or MUSIC
%         zz2{ir,izz} = (R+Ry)/2; % improved MMUSIC cov matrix
%         zz3{ir,izz} = [R Ry]; % improved IMMUSIC cov matrix

        % calculate signal along boresight
        a = ones(antN,1);
        db_psd(ir,izz) = 10*log10(abs(a'*R*a/(a'*a))/df);
        
    end
end