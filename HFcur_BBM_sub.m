%% FFT windows
switch window_FFT
    case 1
        % Blackman Harris window for segment FFTs
        w_BF = window(@blackmanharris,n_BF);
        w_BS = window(@blackmanharris,n_BS);
        w_MU = window(@blackmanharris,n_MU);
        % add other windows here if wanted e.g., case 2...
end

%% Beamform setup
switch window_BF
    case 1
        w_bform = window(@hamming,antN);         % Hamming window
    case 2
        w_bform = ultra(antN, 0.4,1.091,'xmu');  % Ultraspherical window
end

%% Beamscan setup
switch window_BS
    case 1
        w_bscan = window(@hamming,antN);         % Hamming window
    case 2
        w_bscan = ultra(antN, 0.4,1.091,'xmu');  % Ultraspherical window
end

%% directory setup
if ~exist(savedir,'dir') % make save directory if it doesn't exist
    mkdir(savedir)
end

%% start reading sort files
f = dir(fullfile(sortdir,'*.SORT'));
if isempty(f) % check files exist
    error('data directory error')
end

Nfiles = length(f);
for i = 1:Nfiles % for each file
    if Nfiles == 1
        fn = f.name;
    else
        fn = f(i).name; % get filename
    end
    % fn1 = fn(1:11); % YYYYDDDHHMM
    % dirnum = str2double(fn1);

    % save filenames
    fn_short    = fn(1:end-5); % cut off .SORT
    fn_beamIQ   = fullfile(savedir,[fn_short '_beamIQ.mat']);
    fn_beamform = fullfile(savedir,[fn_short '_beamform.mat']);
    fn_beamscan = fullfile(savedir,[fn_short '_beamscan.mat']);
    fn_music    = fullfile(savedir,[fn_short '_music.mat']);

    if strcmp(fn(end-3:end),'SORT')  % only process SORT files

        % read SORT file
        fullname = [sortdir '/' fn]; % full path
        disp(fullname) % display filename

        try
            [WERA,t,R,I,Q]=read_WERA_sort(fullname); % load SORT file
            if sum(isnan(I(:))) > 0 || sum(isnan(Q(:))) > 0
                disp('this file has errors')
                continue
            end
        catch
            disp('this file failed to be read')
            continue
        end

        timer1 = tic; % start timer

        % WERA frequency and range offset
        [r,f_offset] = wera_offsets(R,WERA,ddsclk);

        chirpN      = length(t);
        f_radar     = WERA.FREQ;
        f_bragg     = WERA.fbragg;
        fs          = 1/WERA.RATE;


        % deep water approx, change to full dispersion if you have depth
        % w^2 = g*k
        % lambda = g/(2*pi*f^2);
        %         lambda_bragg = 9.81/2/pi/f_bragg^2;
        lambda_bragg = 299.8/f_radar/2;

        % antenna pattern for ideal monopoles (isotropic) within beam_lim
        [a,theta1,theta_radar] = ant_phases_iso(beam_lim,f_radar,ant_lons,ant_lats,dtheta_beam,r,ant_plot,site,theta_guess);
        % theta = theta1 + theta_radar

        %         % only need antenna pattern within the beam limits 'beam_lim'
        %         itheta      = theta1 >= -beam_lim & theta1 <= beam_lim;
        %         theta111    = theta1(itheta);
        %         theta_grid  = theta(itheta);
        %         zph         = zphi(:,itheta,:);

        %% Beamforming
        if run_beamforming == 1
            Beamform;
        end 

        %% Beamscan 
        if run_beamscan == 1
            Beamscan;
        end

        %% MUSIC
        if run_music == 1
            MUSIC;
        end 

        %% compare Beamforming, Beamscan and MUSIC results from this range cell
        if plot_DOA_analysis == 1
            for ir = 1:length(r)
            % ir = 1;
                sz = 50;
                figure
                clear l1 l2 l3
                % Beamscan DOAs
                l1 = plot(beamscan_data(ir).theta, beamscan_data(ir).vdop,'bs');
                hold on
                plot(beamscan_data(ir).thetan, beamscan_data(ir).vdopn,'bs')
    
                % MUSIC DOAs
                [~,i_snr] = sort(music_data(ir).signals);
                if size(i_snr,1) > 1
                    l2 = scatter(reshape(music_data(ir).thetas(i_snr),[],1),reshape(music_data(ir).vdop(i_snr),[],1),sz,reshape(10*log10(music_data(ir).signals(i_snr)),[],1),'r>');
                else
                    l2 = scatter(reshape(music_data(ir).thetas(i_snr),[],1),repmat(music_data(ir).vdop,size(i_snr,2),1),sz,reshape(10*log10(music_data(ir).signals(i_snr)),[],1),'r>');
                end
                hold on
                [~,i_snr] = sort(music_data(ir).signalsn);
                if size(i_snr,1) > 1
                    scatter(reshape(music_data(ir).thetasn(i_snr),[],1),reshape(music_data(ir).vdopn(i_snr),[],1),sz,reshape(10*log10(music_data(ir).signalsn(i_snr)),[],1),'r>')
                else
                    scatter(reshape(music_data(ir).thetasn(i_snr),[],1),repmat(music_data(ir).vdopn,size(i_snr,2),1),sz,reshape(10*log10(music_data(ir).signalsn(i_snr)),[],1),'r>')
                end
       
             
                % Beamforming reults
                l3 = plot(theta1,V(ir,:),'ks');
    
                title(['range cell ' num2str(ir)])
                xlabel('\theta from boresight')
                ylabel('current (m/s)')
                
                if ~isempty(l1) && ~isempty(l2) && ~isempty(l3)
                    legend([l1 l2 l3],'Beamscan','MUSIC','Beamform')
                end
                drawnow
            end
        end

        %% with colored signal strength
        if plot_DOA_analysis == 2 %  SNR for Beamscan and signal for MUSIC
            for ir = 1:length(r)
            % ir = 1;
                sz = 50;
                figure
                
                clear l1 l2 l3
                % Beamscan DOAs
                ax1 = axes;
                [~,i_snr] = sort(beamscan_data(ir).snr);
                l1 = scatter(beamscan_data(ir).theta(i_snr),beamscan_data(ir).vdop(i_snr),sz,beamscan_data(ir).snr(i_snr),'filled','o');
                hold on
                [~,i_snr] = sort(beamscan_data(ir).snrn);
                scatter(beamscan_data(ir).thetan(i_snr),beamscan_data(ir).vdopn(i_snr),sz,beamscan_data(ir).snrn(i_snr),'filled','o');
     
                title(['range cell ' num2str(ir)])
                xlabel('\theta from boresight')
                ylabel('current (m/s)')
                
                % MUSIC DOAs
                ax2 = axes;
                [~,i_snr] = sort(music_data(ir).signals);
                if size(i_snr,1) > 1
                    l2 = scatter(reshape(music_data(ir).thetas(i_snr),[],1),reshape(music_data(ir).vdop(i_snr),[],1),sz,reshape(10*log10(music_data(ir).signals(i_snr)),[],1),'filled','>');
                else
                    l2 = scatter(reshape(music_data(ir).thetas(i_snr),[],1),repmat(music_data(ir).vdop,size(i_snr,2),1),sz,reshape(10*log10(music_data(ir).signals(i_snr)),[],1),'filled','>');
                end
                hold on
                [~,i_snr] = sort(music_data(ir).signalsn);
                if size(i_snr,1) > 1
                    scatter(reshape(music_data(ir).thetasn(i_snr),[],1),reshape(music_data(ir).vdopn(i_snr),[],1),sz,reshape(10*log10(music_data(ir).signalsn(i_snr)),[],1),'filled','>')
                else
                    scatter(reshape(music_data(ir).thetasn(i_snr),[],1),repmat(music_data(ir).vdopn,size(i_snr,2),1),sz,reshape(10*log10(music_data(ir).signalsn(i_snr)),[],1),'filled','>')
                end
                
                % Beamforming reults
                l3 = plot(theta1,V(ir,:),'ks');
                
                %Link them together
                linkaxes([ax1,ax2]);
                %Hide the top axes
                ax2.Visible = 'off';
                ax2.XTick = [];
                ax2.YTick = [];
                %Give each one its own colormap
                colormap(ax1,'parula')
                colormap(ax2,'parula')
                
                if ~isempty(l1) && ~isempty(l2) && ~isempty(l3)
                    legend([l1 l2 l3],'Beamscan','MUSIC','Beamform')
                end
                drawnow
                pause(0.25)
            end
        end

    end % for SORT file
end % for each file













