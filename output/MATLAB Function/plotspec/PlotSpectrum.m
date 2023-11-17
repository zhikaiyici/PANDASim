function [counts, edges] = PlotSpectrum(spectrumData,channalWidth)
% channalwidth MeV

h = histogram(spectrumData, 'BinWidth', channalWidth,...
    'DisplayStyle', 'stairs', 'EdgeColor', 'k');
counts = h.Values;
edges = h.BinEdges;
% ax = gca;
set(gca, 'fontname', 'times new roman', 'xgrid', 'on', 'ygrid', 'on');
maxCounts = max(counts);
set(gca, 'ylim', [0, 1.1 * maxCounts]);
% set(gca, 'yscale', 'log');
% minedep = min(spectrumdata);
maxEdep = max(spectrumData);
set(gca, 'xlim', [0, 1.1 * maxEdep]);
xlabel("Energy (MeV) / Photoelectron");
ylabel("Counts");

end
