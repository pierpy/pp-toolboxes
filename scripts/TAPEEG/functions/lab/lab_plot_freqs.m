function F = lab_plot_freqs(SpectAllM,SpectAllF,IFREQ,Filename)

if ~exist('IFREQ','var')
    IFREQ = [];
end
Idx50 = find(SpectAllF>=50,1,'first');
SpectAllM = SpectAllM(:,1:Idx50);
SpectAllF = SpectAllF(1:Idx50);

Spect = median(SpectAllM,1);
Freqlabel = cell(1,length(SpectAllF));
for i = 1:length(SpectAllF)
    if round(SpectAllF(i)) == SpectAllF(i)
        Freqlabel(i) = num2cell(SpectAllF(i));
    end
end
tmp = find(~cellfun(@isempty,Freqlabel));
if ~isempty(tmp) & length(tmp) >25
    for i = 2:2:length(tmp)
        Freqlabel{1,tmp(i)} = [];
    end
elseif isempty(tmp) | length(tmp) < 5
    Freqlabel = cell(1,length(SpectAllF));
    for i = 1:floor(length(SpectAllF)/25):length(SpectAllF)
        Freqlabel{1,i} = SpectAllF(i);
    end
end
clearvars tmp

if exist('Filename','var') & ~isempty(Filename)
    F = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Visible','off');
else
    F = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off');
end
area(Spect);
set(gca,'XTick',1:length(Freqlabel));
set(gca,'XTickLabel',Freqlabel,'FontName','Times','fontsize',9);

if isfield(IFREQ,'Bands') & ~isempty(IFREQ.Bands)
    Lsize = max(Spect);
    for i = 1:size(IFREQ.Bands,1)
        hold on
        Idx = find(SpectAllF >= IFREQ.Bands{i,2},1);
        Idx2 = find(SpectAllF >= IFREQ.Bands{i,3},1);
        switch IFREQ.Bands{i,1}
            case 'Delta'
                patch([Idx Idx2 Idx2 Idx],[0 0 Lsize Lsize],'b','FaceAlpha',0.3);
            case 'Theta'
                patch([Idx Idx2 Idx2 Idx],[0 0 Lsize Lsize],'y','FaceAlpha',0.3);
            case 'Alpha1'
                patch([Idx Idx2 Idx2 Idx],[0 0 Lsize Lsize],'r','FaceAlpha',0.3);
            case 'Alpha2'
                patch([Idx Idx2 Idx2 Idx],[0 0 Lsize Lsize],'r','FaceAlpha',0.3);
            case 'Beta'
                patch([Idx Idx2 Idx2 Idx],[0 0 Lsize Lsize],'g','FaceAlpha',0.3);    
            case 'Alpha1a'
                patch([Idx Idx2 Idx2 Idx],[0 0 0.6*Lsize 0.6*Lsize],'r','FaceAlpha',0.3);    
            case 'Alpha1b'
                patch([Idx Idx2 Idx2 Idx],[0 0 0.6*Lsize 0.6*Lsize],'r','FaceAlpha',0.3);
            case 'Beta1'
                patch([Idx Idx2 Idx2 Idx],[0 0 0.6*Lsize 0.6*Lsize],'g','FaceAlpha',0.3);    
            case 'Beta2'
                patch([Idx Idx2 Idx2 Idx],[0 0 0.6*Lsize 0.6*Lsize],'g','FaceAlpha',0.3);     
        end
    end
end

if exist('Filename','var') & ~isempty(Filename)
    lab_print_figure(Filename,F);
end