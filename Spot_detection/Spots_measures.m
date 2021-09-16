
close all
clc

% NOTE: this script requires the following measures to have been made and exported
% ObjectInclusionCount -> Containing: Bacteria; To count: Spot_detection -> Name: SpotCount
% SpineCoordinates -> Set SpineLength to parent: false
% SpineFeatures -> Bacteria; Scaled: pixel
% Optional Measurement:
% ObjectFeatures -> Object class: Bacteria -> Feature: Mean -> Intensity: SOS_signal -> Name: MeanSOS


Bacmman_folder = '/media/daniel/HDD Daniel/Daniel Thédié/BACMMAN/Timelapse/'; % Bacmman working directory
dataset_name = '210914_DT7'; % Bacmman dataset name
prefix = 'Im_'; % Prefix found in all images before the number (e.g. for Im1, Im2,... set prefix = 'Im';)
files_folder = '/media/daniel/HDD Daniel/Daniel Thédié/Timelapse/210914/DT7/'; % Folder where the original images are stored (only used to retrieve their timestamp)
bf_keyword = '_w10 Brightfield'; % Keyword to identify brightfield images

heatmap_video = 0; % Set to 1 to make a video of the heatmaps of spot positions in time
fps = 20; % Number of frames per second for the video

%% End of input

% Import data
dataCells = readtable([Bacmman_folder dataset_name '/' dataset_name '_0.csv'], 'TreatAsEmpty', 'NA');
dataSpots = readtable([Bacmman_folder dataset_name '/' dataset_name '_1.csv'], 'TreatAsEmpty', 'NA');

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

% Fetch timestamps for each FOV
cd(files_folder)
listing = dir(['*' bf_keyword '*']);
uCellIdx = unique(dataCells.TrueIdx);
for i = 1:length(uCellIdx)
    bfInfo = imfinfo(listing(i).name);
    dataCells.Timestamp(dataCells.TrueIdx == uCellIdx(i)) = datetime(bfInfo(1).DateTime, 'InputFormat', 'yyyyMMdd HH:mm:ss.SSS');
end
dataCells.Time = dataCells.Timestamp - dataCells.Timestamp(1);
timepoints = unique(dataCells.Time);

% Measures
fovNum = unique(dataCells.TrueIdx);
nCellsFOV = grpstats(dataCells.Idx, dataCells.TrueIdx, 'max') +1;
nCells = height(dataCells);
spotFrac = sum(dataCells.SpotCount > 0)/nCells;
spotFracFOV = grpstats(dataCells.SpotCount > 0, dataCells.TrueIdx, 'sum')./nCellsFOV;
zeroSpot = grpstats(dataCells.SpotCount == 0, dataCells.TrueIdx, 'sum')./nCellsFOV;
oneSpot = grpstats(dataCells.SpotCount == 1, dataCells.TrueIdx, 'sum')./nCellsFOV;
twoPlusSpot = grpstats(dataCells.SpotCount >= 2, dataCells.TrueIdx, 'sum')./nCellsFOV;
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    sosFOV = grpstats(dataCells.MeanSOS, dataCells.TrueIdx, 'mean');
end
cellLenFOV = grpstats(dataCells.SpineLength, dataCells.TrueIdx, 'mean');


% Normalise spot positions in the cell
% SpineCurvilinearCoord: [0:cellLength]
% SpineRadialCoord: distance from spine (pixels)

dataSpots.normXpos = dataSpots.SpineCurvilinearCoord./dataSpots.SpineLength -0.5; % Normalise to cell length & centre
dataSpots.normYpos = dataSpots.SpineRadialCoord./dataSpots.SpineRadius; % Normalise to cell radius





%% Figures

fprintf('Average number of cells per FOV: %.0f +/- %.0f\n', mean(nCellsFOV), std(nCellsFOV))
fprintf('Proportion of cells with spots (all FOVs): %.0f%% (n=%.0f)\n', spotFrac*100, nCells)
fprintf('Average cell length: %.0f +/- %.0f pixels\n', mean(dataCells.SpineLength), std(dataCells.SpineLength))
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    fprintf('Average SOS signal: %.0f +/- %.0f\n', mean(dataCells.MeanSOS), std(dataCells.MeanSOS))
end


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

% Cell length vs time
length_movmean = movmean(cellLenFOV, 50);
figure('Color', 'white')
hold on
scatter(minutes(timepoints), cellLenFOV, 15, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerFaceColor', 'b');
plot(minutes(timepoints), length_movmean, 'LineWidth', 1.5, 'Color', 'r');
hold off
box on
xlabel('Time (min)')
ylabel('Average cell length (pixels)')



% Histogram of SOS signal
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    figure('Color', 'white')
    histogram(dataCells.MeanSOS, 'Normalization', 'probability', 'BinWidth', 200)
    xlabel('Mean SOS signal per cell')
    ylabel('PDF')
else
    fprintf('No SOS data found.\n')
end

% SOS signal vs time
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    sos_movmean = movmean(sosFOV, 50);
    figure('Color', 'white')
    hold on
    scatter(minutes(timepoints), sosFOV, 15, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerFaceColor', 'b');
    plot(minutes(timepoints), sos_movmean, 'LineWidth', 1.5, 'Color', 'r');
    hold off
    box on
    xlabel('Time (min)')
    ylabel('Average SOS signal')
end


% Cell length/SOS signal correlation
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    figure('Color', 'white')
    scatter(cellLenFOV, sosFOV, 15, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', 'b');
    box on
    xlabel('Cell length (pixels)')
    ylabel('SOS signal (AU)')
end


% SOS signal vs time
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    figure('Color', 'white')
    hold on
    y = zeros(3, length(timepoints));
    for i = 1:length(timepoints)
        y(1, i) = mean(dataCells.MeanSOS(dataCells.TrueIdx == i & dataCells.SpotCount == 0));
        y(2, i) = mean(dataCells.MeanSOS(dataCells.TrueIdx == i & dataCells.SpotCount == 1));
        y(3, i) = mean(dataCells.MeanSOS(dataCells.TrueIdx == i & dataCells.SpotCount >= 2));
    end
    sosMovMeanZero = movmean(y(1,:), 50, 'omitnan');
    sosMovMeanOne = movmean(y(2,:), 50, 'omitnan');
    sosMovMeanTwo = movmean(y(3,:), 50, 'omitnan');
    scatter(minutes(timepoints), y(1,:), 15, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'k');
    scatter(minutes(timepoints), y(2,:), 15, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'b');
    scatter(minutes(timepoints), y(3,:), 15, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'm');
    plot(minutes(timepoints), sosMovMeanZero, 'LineWidth', 1.5, 'Color', 'k');
    plot(minutes(timepoints), sosMovMeanOne, 'LineWidth', 1.5, 'Color', 'b');
    plot(minutes(timepoints), sosMovMeanTwo, 'LineWidth', 1.5, 'Color', 'm');
    hold off
    box on
    xlabel('Time (min)')
    ylabel('Average SOS signal')
    legend('No spot','One spot','2 or more spots')
end




% Number of spots per cell histogram
figure('Color', 'white')
histogram(dataCells.SpotCount, 'Normalization', 'probability', 'BinWidth', 1,...
    'BinEdges', -0.5:1:max(dataCells.SpotCount)+0.5)
xlabel('Number of spots per cell')
ylabel('PDF')


% Time evolution of the number of spots per cell
zeroSpot_mean = movmean(zeroSpot, 50);
oneSpot_mean = movmean(oneSpot, 50);
twoSpot_mean = movmean(twoPlusSpot, 50);
figure('Color', 'white')
hold on
scatter(minutes(timepoints), zeroSpot, 15, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'k');
scatter(minutes(timepoints), oneSpot, 15, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'b');
scatter(minutes(timepoints), twoPlusSpot, 15, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'm');
plot(minutes(timepoints), zeroSpot_mean, 'LineWidth', 1.5, 'Color', 'k');
plot(minutes(timepoints), oneSpot_mean, 'LineWidth', 1.5, 'Color', 'b');
plot(minutes(timepoints), twoSpot_mean, 'LineWidth', 1.5, 'Color', 'm');
hold off
box on
xlabel('Time (min)')
ylabel('Proportion of cells containing 0, 1 or more spots')
ylim([0 1.3])
legend('No spot','One spot','2 or more spots')




% Histogram of fraction of cells with spots per FOV
figure('Color', 'white')
histogram(spotFracFOV, 'Normalization', 'probability', 'BinWidth', 0.02)
xlabel('Fraction of cells with spots per FOV')
ylabel('PDF')
xlim([0 1])



% Time evolution of the fraction of cells with spots per FOV
spotFrac_mean = movmean(spotFracFOV, 50);
figure('Color', 'white')
hold on
scatter(minutes(timepoints), spotFracFOV, 15, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerFaceColor', 'b');
plot(minutes(timepoints), spotFrac_mean, 'LineWidth', 1.5, 'Color', 'r');
hold off
box on
xlabel('Time (min)')
ylabel('Proportion of cells containing an immobile spot')
ylim([0 1])
legend('Individual FOVs', '50-points moving average')


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


if heatmap_video
    
    if exist('F', 'var')
        clear F
    end
    
    figure('Color', 'white')
    for i = 1:length(timepoints)-20
        
        xypos = [dataSpots.normXpos(dataSpots.TrueIdx >= i & dataSpots.TrueIdx < i+20)...
            dataSpots.normYpos(dataSpots.TrueIdx >= i & dataSpots.TrueIdx < i+20)];
        mapData = xypos((xypos(:,1) < 2 & xypos(:,1) > -2) & (xypos(:,2) < 2 & xypos(:,2) > -2), :);
        
        hist3(mapData,'CdataMode','auto', 'edges', {-1:0.075:1 -2:0.15:2})
        hold on
        text(0.15, 1.5, [num2str(round(minutes(timepoints(i+10)))) ' min'],...
            'Color', 'w', 'FontWeight', 'b', 'FontSize', 25)
        hold off
        colormap(jet)
        title('Position of bright spots (all cells)')
        xlabel('% of major axis')
        ylabel('% of minor axis')
        colorbar
        view(2)
        xlim([-1 1])
        ylim([-1*2 1*2])
        caxis([0 50])
        
        F(i) = getframe(gcf) ;
        drawnow
    end
    
    % create the video writer with specified fps
    writerObj = VideoWriter([Bacmman_folder dataset_name '/Heatmap_video.avi']);
    writerObj.FrameRate = fps;
    writerObj.Quality = 100;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=2:length(F)
        % convert the image to a frame
        frame = F(i) ;
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    
end






































