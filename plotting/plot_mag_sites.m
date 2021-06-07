function plot_mag_sites(b,tidx,is,print_flag)

    subplot(2,1,1);
    plot(b(is).times(tidx),b(is).x(tidx)-nanmean(b(is).x(tidx)),'-b'); hold on
    plot(b(is).times(tidx),b(is).y(tidx)-nanmean(b(is).y(tidx)),'-r');
    ylabel('B (nT)')
    title(['Magnetic Time Series for ',upper(b(is).site)]);

    subplot(2,1,2);
    loglog(b(is).f,abs(b(is).X(1:length(b(is).f))),'-b'); hold on
    loglog(b(is).f,abs(b(is).Y(1:length(b(is).f))),'-r');
    xlabel('Frequency (Hz)'); ylabel('abs(B)')
    title(['Magnetic Spectrum for ',upper(b(is).site)]);

    if print_flag == 1
        print_figure('',['Mag_Site_',b(is).site]);
    end