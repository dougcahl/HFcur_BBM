function [theta00,thetas00,signal1,signals1] = music_single_bragg_peakc1_sub(ant_pattern_phi,theta1,nmax,peak_thresh,Z,plot_MUSIC_DOAs,tname)

theta00 = nan;
thetas00 = nan(1,nmax);
signal1 = nan;
signals1 = nan(1,nmax);

tsize = size(Z);
if tsize(1) == tsize(2)
    [V,D] = eig(Z); % diagonal matrix D of eigenvalues and V columns eigenvectors
    D = diag(D);
else
    [V,D,~] = svd(Z);
    D = diag(D);
    D = D(1:tsize(1));
end

[D] = sort(D);

for n_dir = nmax:-1:1
    EN = V(:,1:end-n_dir); % leaving n_dir biggest eigenvectors
    ES = V(:,end-n_dir+1:end); % signal eig vectors, n_dir biggest eigenvectors
    % 1/d^2 = Y* EN EN* Y
    % 1/( a*(t) En EN* a(t) )
    pt = nan(size(theta1));
    for ien = 1:length(theta1)
        at = ant_pattern_phi(:,ien); % a(theta)
        d2 = at' * (EN * EN') * at;
        pt(ien) = 1/abs(d2);
    end
    % get peak max
    [ptm,pti] = max(pt);
    pt = pt(:);

    if plot_MUSIC_DOAs == 1
        hold off
        subplot(511)
        hold off
        semilogy(D,'bs')
        xlim([0 13])
        hold on
        title([tname ' eig vals'])
    
        subplot(5,1,n_dir+1)
        hold off
        plot(theta1,pt)
        hold on
        plot(theta1(pti),ptm,'rs')%,'markerfacecolor','r')
        title([num2str(n_dir) ' DOAs'])
    end
    
    
    [pks,dzdt]=pksfinder1(10*log10(pt),peak_thresh);
    ipks=find(diff(dzdt)==-2)+1;
    pks_cnt = length(ipks);
    apks = pt(ipks);
    alocs = ipks;
    if plot_MUSIC_DOAs == 1
        % draw peaks
        plot(theta1(alocs),apks,'v')
    end
    
    if pks_cnt == n_dir || n_dir == 1 % save multiple peaks here
        
        apks = pt(ipks);
        alocs = ipks;
        [apks,alocsi] = sort(apks,'descend');
        alocs = alocs(alocsi);
        for ialocs = 1:n_dir
            if isempty(alocs)
                thetas00(ialocs) = nan;
                signals1(ialocs) = nan;
                break
            end
            thetas00(ialocs) = theta1(alocs(ialocs));
            at1 = ant_pattern_phi(:,alocs(ialocs)); % a(theta)
            G = at1'*ES;
            eig_matrix = D(end-n_dir+1:end);
            eig_matrix = diag(eig_matrix);
            
            signals1_matrix = real(G'\eig_matrix/G);
            signals1_snrs = diag(signals1_matrix);
            signals1(ialocs) = signals1_snrs; 
        end
        break % break here (save unnecessary processing)
    end
end % for n_dir

if plot_MUSIC_DOAs == 1
    drawnow
end

% removed calculations for single DOA but you can replace them down here if
% you really want a single DOA calculation for every frequency (not
% accurate). It's the same code as above just for n_dir=1 and instead of
% running using the pksfinder routine you simply choose the location of the
% maximum value of the pseudo-spectrum for 1 DOA