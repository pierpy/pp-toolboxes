% Show results of spectral analysis (used by lab_load_eeg)
%
% written by F. Hatz

function lab_show_spectraresults(Spectra,Mappings)

mode = questdlg('Mean or Median','Show spectral results','Mean','Median','Median');

freqlabel = cell(1,length(Spectra.Power_freqs));
for i = 1:length(Spectra.Power_freqs)
    if round(Spectra.Power_freqs(i)) == Spectra.Power_freqs(i)
        freqlabel(i) = num2cell(Spectra.Power_freqs(i));
    end
end
if strcmp(mode,'Median')
    plot_spectra(median(Spectra.Power_median,1),freqlabel,1,'Median Spectral density');
else
    plot_spectra(mean(Spectra.Power_mean,1),freqlabel,1,'Mean Spectral density');
end

if isfield(Spectra,'BandPower_freqs')
    freqlabel = Spectra.BandPower_freqs;
    if strcmp(mode,'Median')
        plot_spectra(Spectra.BandPower_median',freqlabel,2,'BandPower');
    else
        plot_spectra(Spectra.BandPower_mean',freqlabel,2,'BandPower');
    end
end
if isfield(Spectra,'BP_Mapping_mean')
    freqlabel = Spectra.BandPower_freqs;
    if exist('Mappings','var') & isfield(Mappings,'mappingstitle')
        legends = Mappings.mappingstitle;
    else
        legends = [];
    end
    if strcmp(mode,'Median')
        plot_spectra(Spectra.BP_Mapping_median',freqlabel,2,'BandPower-Mappings',legends);
    else
        plot_spectra(Spectra.BP_Mapping_mean',freqlabel,2,'BandPower-Mappings',legends);
    end
end

end

function plot_spectra(SpectAll,freqlabel,mode,Title,legends)

figure('Color',[1 1 1]);
if ~exist('mode','var') | mode == 1
    area(SpectAll);
else
    plot(SpectAll);
    if exist('legends','var') & ~isempty(legends)
        legend(legends);
    end
end
set(gca,'XTick',1:length(freqlabel));
set(gca,'XTickLabel',freqlabel,'FontName','Times','fontsize',9);
if exist('Title','var')
    title(Title);
    set(gcf,'NumberTitle','off','Name',Title);
end

end