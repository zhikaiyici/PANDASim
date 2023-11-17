
name = ls;
fileName = 'moduleEdepDelayGd(3x3_1e+06_neutrino_100%_5).txt';%name(18,:);
% fileName = 'moduleEdepPrompt(3x3_1e+05_Co60_CENTER).txt';%name(18,:);
arraySize = 3;
channalWidth = 0.05; % channalwidth MeV / pheNum

spectraData = load(fileName);
spectraData = ReshapeDataMatrix(arraySize, spectraData);
spectraData(spectraData < channalWidth) = 0;

% % spectraData = spectraData .* 1e7;
% sigma = sqrt(spectraData .* 9 ./ 8);
% for ii = 1:numel(spectraData)
%     spectraData(ii) = normrnd(spectraData(ii),sigma(ii));
% end
% spectraData(spectraData < 0) = 0;

logic = spectraData >= channalWidth;
if arraySize > 1
    triggerEvents = sum(sum(logic));
else
    triggerEvents = logic;
end
figure('Name', ['Trigger_', fileName]);
h = histogram(triggerEvents, 'BinEdges', 0.5:arraySize * arraySize + 0.5, 'Visible', 'off');
triggerEvents = h.Values;
triggerNum = h.BinEdges - 0.5;
triggerNum(triggerNum == 0) = []; 
bar(triggerNum, triggerEvents);

if arraySize > 1
    totalEdep = sum(sum(spectraData));
else
    totalEdep = spectraData;
end
totalEdep(totalEdep == 0) =[];
figure('Name', ['Total_', fileName])
PlotSpectrum(totalEdep, channalWidth);

figure('Name', fileName)
subNum = 0;
for ii = 1:arraySize
    for jj = 1:arraySize
        moduleEdep = spectraData(ii, jj, :);
        moduleEdep(moduleEdep == 0) = [];
        subNum = subNum + 1;
        subplot(arraySize,arraySize, subNum);
        PlotSpectrum(moduleEdep, channalWidth);
    end
end
