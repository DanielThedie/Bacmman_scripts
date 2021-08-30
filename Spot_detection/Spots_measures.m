
close all
clc

% NOTE: this script requires the following measures to have been made and exported
% ObjectInclusionCount -> Containing: Bacteria; To count: Spot_detection -> Name: SpotCount
% SpineCoordinates -> Set SpineLength to parent: false
% SpineFeatures -> Bacteria; Scaled: pixel
% Facultative Measurement:
% ObjectFeatures -> Object class: Bacteria -> Feature: Mean -> Intensity: SOS_signal -> Name: MeanSOS


Bacmman_folder = '/media/daniel/HDD Daniel/Daniel Thédié/BACMMAN/'; % Bacmman working directory
dataset_name = '210827_DT7'; % Bacmman dataset name
prefix = 'Im'; % Prefix found in all images before the number (e.g. for Im1, Im2,... set prefix = 'Im';)


%% End of input

% Import data
dataCells = readtable([Bacmman_folder dataset_name '/' dataset_name '_0.csv']);
dataSpots = readtable([Bacmman_folder dataset_name '/' dataset_name '_1.csv']);

% Sort datasets chronologically
for i = 1:height(dataCells)
    name = dataCells.Position(i);
    if ~isempty(prefix)
        name = char(name);
        idx = name(strfind(name, prefix)+length(prefix):end);
        dataCells.TrueIdx(i) = str2double(idx);
    else
        dataCells.TrueIdx(i) = name;
    end
end
dataCells = sortrows(dataCells, 'TrueIdx', 'ascend');
for i = 1:height(dataSpots)
    name = dataSpots.Position(i);
    if ~isempty(prefix)
        name = char(name);
        idx = name(strfind(name, prefix)+length(prefix):end);
        dataSpots.TrueIdx(i) = str2double(idx);
    else
        dataSpots.TrueIdx(i) = name;
    end
end
dataSpots = sortrows(dataSpots, 'TrueIdx', 'ascend');


% Measures
fovNum = unique(dataCells.TrueIdx);
nCellsFOV = grpstats(dataCells.Idx, dataCells.TrueIdx, 'max');
nCells = height(dataCells);
spotFrac = sum(dataCells.SpotCount > 0)/nCells;
spotFracFOV = grpstats(dataCells.SpotCount > 0, dataCells.TrueIdx, 'sum')./nCellsFOV;


% Normalise spot positions in the cell
% SpineCurvilinearCoord: [0:cellLength]
% SpineRadialCoord: distance from spine (pixels)

dataSpots.normXpos = dataSpots.SpineCurvilinearCoord./dataSpots.SpineLength -0.5; % Normalise to cell length & centre
dataSpots.normYpos = dataSpots.SpineRadialCoord./dataSpots.SpineRadius; % Normalise to cell radius





%% Figures

fprintf('Average number of cells per FOV: %.0f +/- %.0f\n', mean(nCellsFOV), std(nCellsFOV));
fprintf('Proportion of cells with spots (all FOVs): %.0f%% (n=%.0f)\n', spotFrac*100, nCells)


% Number of cells recorded per FOV
figure('color','white')
scatter(fovNum, nCellsFOV, 15, 'filled')
box on
xlabel('FOV number')
ylabel('Number of cells recorded')
ylim([0 1.2*max(nCellsFOV)])


% Histogram of cell lengths
figure('Color', 'white')
histogram(dataCells.SpineLength, 'Normalization', 'probability')
xlabel('Cell length (pixels)')
ylabel('PDF')


% Histogram of SOS signal
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    figure('Color', 'white')
    histogram(dataCells.MeanSOS, 'Normalization', 'probability')
    xlabel('Mean SOS signal per cell')
    ylabel('PDF')
else
    fprintf('No SOS data found.\n')
end

% Number of spots per cell histogram
figure('Color', 'white')
histogram(dataCells.SpotCount, 'Normalization', 'probability', 'BinWidth', 1,...
    'BinEdges', -0.5:1:max(dataCells.SpotCount)+0.5)
xlabel('Number of spots per cell')
ylabel('PDF')



% Histogram of fraction of cells with spots per FOV
figure('Color', 'white')
histogram(spotFracFOV, 'Normalization', 'probability', 'BinWidth', 0.02)
xlabel('Fraction of cells with spots per FOV')
ylabel('PDF')
xlim([0 1])



% Time evolution of the fraction of cells with spots per FOV
figure('Color', 'white')
hold on
for i = 1:length(nCellsFOV)
    s1 = scatter(i, spotFracFOV(i), 15, 'filled');
    s1.MarkerFaceAlpha = nCellsFOV(i)/max(nCellsFOV);
    s1.MarkerFaceColor = 'b';
end
hold off
box on
xlabel('FOV number')
ylabel('Proportion of cells containing an immobile spot')
ylim([0 1])



% Heatmap of spot positions in the cell (normalised)
xypos = [dataSpots.normXpos dataSpots.normYpos];
mapData = xypos((xypos(:,1) < 2 & xypos(:,1) > -2) & (xypos(:,2) < 2 & xypos(:,2) > -2), :);

figure('Color', 'white')
hist3(mapData,'CdataMode','auto', 'edges', {-1:0.075:1 -2:0.15:2})
colormap(jet)
title('Position of bright spots (all cells)')
xlabel('% of major axis')
ylabel('% of minor axis')
colorbar
view(2)
xlim([-1 1])
ylim([-1*2 1*2])









































