function [freq,db_psd,zz] = fft_cov(I,Q,w,fs)
%#codegen
% [db,zz,zzmat] = function fft_cov(data,wn3,n_seg,p_overlap)
%
antN = size(I,1);
rN = size(I,2);
chirpN = size(I,3);
u = 1/length(w)*sum(w.^2);
freq = 0;

    % FFT calcs
    db_psd = nan(antN,rN,chirpN);
    zz = complex(db_psd);
%     zz_mat = cell(rN,chirpN+1);
    for ir = 1:rN
        for in = 1:antN
            %%%%%%%%%%%%%%%%%%%%%%%%%%% full length FFT
            data = squeeze(I(in,ir,:)) + 1i*squeeze(Q(in,ir,:));
            data = data.*w;
            [freq,~,~,psd,~,zz1] = fft_doug_phase_2(data, fs);
            %             energy      = (1/u)*energy;
            psd = (1/u)*psd;
            %             amp = sqrt(1/u)*amp;
            zz1 = sqrt(1/u)*zz1;
            db_psd1 = 10*log10(psd);
            db_psd(in,ir,:) = db_psd1;
            zz(in,ir,:) = zz1;
        end
        % now we can make matrices using all 12 antennas for full length FFT
%         for izz = 1:chirpN+1
%             zz11 = zz(:,ir,izz);
%             zz12 =  zz11*zz11';
%             zz_mat{ir,izz} = zz12;
%         end
    end

