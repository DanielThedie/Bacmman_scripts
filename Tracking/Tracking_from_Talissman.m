
close all
clc

addpath '../Util'
addpath 'Tracking_util' % These two do not need to be updated if you've downloaded the whole Github repository

% Dependencies
% Get_cells.m, findpeak.m, track.m

% Folder containing streams
folder = '/media/daniel/HDD Daniel/Daniel Thédié/Tracking/210820/AL103_100pc/';
streamName = 'Stream'; % Name of "stream" files to be processed

% Experimental parameters
param.Xtime = 150;     % objective magnification
param.pixel_size = 16; % Chip pixel size

% Segmentation file info
segFile_param.folder = '/media/daniel/HDD Daniel/Daniel Thédié/BACMMAN/Tracking/210820_AL103_100pc/';
segFile_param.fileName = 'Cells.h5';
segFile_param.cellsFeature = 'Cells';

% Peak detection
param.thrfpeak = 10;   % Threshold for peak detection
param.pnoise = 1;      % noise in peakfind function
param.psize = 4;       % size group of pixel in peakfind function
param.pgauss = 5;      % size for gaussian fit peak size in centfind

% Tracking
maxdisp = 5;           % Maximum displacement
minTrackLength = 4;    % Minimum length of tracks (shorter tracks will be removed)

para.mem = 0; % Memory parameter; track won't be cut if the molecule disappears for "mem" frames
para.good = 2;
para.dim = 2;
para.quiet = 1;


%% End of input

trackData = [];
skipped = 0;

cd(folder)
fovs = dir([streamName '*']);
% Sort FOVs so they're in the correct order
tmp = struct2table(fovs);
[~, idx] = sortrows(tmp, 'datenum', 'ascend');
fovs = fovs(idx);

maxParticleID = 0;
maxCellID = 0;
cells_numb = zeros(length(fovs), 1);

for i = 1:length(fovs) % Loop on all FOVs
    
    fprintf('Processing FOV %.0f/%.0f\n', i, length(fovs));
    
    % Import stream
    imgfile = [folder fovs(i).name];
    iminfo = imfinfo(imgfile);
    nx = iminfo.Height;
    ny = iminfo.Width;
    nz = length(iminfo);
    zim = zeros(nx, ny, nz);
    for j = 1:nz
        zim(:,:,j) = imread(imgfile,'Index',j);
    end
    
    % Peak detection
    
    % Columns in the peaks variable:
    % 1,2: x and y width of gaussian fit
    % 3: rotation angle
    % 4: offset
    % 5,6: x and y coordinates
    % 7: peak intensity
    % 8: frame number
    
    param.nstacks = nz;
    peaks = findpeak(zim,param);
    
    % Inversion of the X and Y axis to match the coordinate system of the segmentation
    tmp = peaks(:, 5);
    peaks(:,5) = peaks(:,6);
    peaks(:,6) = tmp;
    
    % Get segmentation info (hdf5 file from Talissman)
    segmentation = Get_cells(segFile_param, i);
    props = {'Area','Centroid','MajorAxisLength','MinorAxisLength', 'PixelList'};
    cells = regionprops(segmentation,props);
    cells_numb(i) = length(cells);
    
    clc
    fprintf('%.0f cells on this FOV\n', length(cells))
    
    
    for j = 1:length(cells) % Loop over all cells
        
        fprintf('%.0f/%.0f\n', j, length(cells));
        
        % Select peaks corresponding to the current cell
        peakSel = peaks(ismember(round(peaks(:,[5,6])), cells(j).PixelList, 'row'), :);
        
        % Tracking
        
        % Columns in trackData:
        % 1, 2: x and y coordinates (pixels, converted to nm straight after tracking)
        % 3: frame number
        % 4: particle ID
        % 5: peak intensity (added after the tracking routine)
        % 6: FOV number
        % 7: cell ID
        
        skip = 0;
        try
            tracks = track(peakSel(:,[5,6,8]), maxdisp, para);
        catch
            [~, outliers] = rmoutliers(peakSel(:,8));
            peakSel = peakSel(~outliers,:);
            try
                tracks = track(peakSel(:,[5,6,8]), maxdisp, para);
            catch
                skip = 1;
                skipped = skipped +1;
                disp('Skipping cell')
            end
        end
        
        if ~skip
            
            disp('Tracking complete')
            % Add peak intensity variable (from peaks) to tracks
            [~, iind] = ismember(tracks(:,1:2),peakSel(:,5:6), 'row');
            tracks(:,5) = peaks(iind,7);
            
            % Convert pixels to nanometres
            tracks(:,1:2) = tracks(:,1:2) * param.pixel_size/param.Xtime *1000;
            
            % Remove short tracks (threshold defined in input)
            d = [true, diff(tracks(:,4).') ~= 0, true];  % TRUE if values change
            n = diff(find(d));                           % Number of repetitions
            Y = repelem(n, n);
            tracks = tracks(Y > minTrackLength,:);
            
            % Increment particle ID so that different cells don't have the same particle ID
            if ~isempty(tracks)
                tracks(:,4) = tracks(:,4) + maxParticleID;
                maxParticleID = max(tracks(:,4));
            end
            
            % Include FOV number
            tracks(:,6) = i*ones(length(tracks), 1);
            
            % Include cell number and increment similarly to particle ID
            tracks(:,7) = j + maxCellID;
            
            fprintf('Removed %.0f/%.0f tracks that were shorter than %.0f frames\n', sum(Y < minTrackLength), length(Y), minTrackLength)
            
            % Update trackData
            trackData = [trackData; tracks];
            
        end
        
    end
    
    maxCellID = max(trackData(:,7));
    
    fovTracks = trackData(trackData(:,6) == i,:);
    
    % List of track numbers for current FOV
    track_nums = unique(fovTracks(:,4));
    
    % Show tracks on maximum projection
    xy_pix = [fovTracks(:,2), fovTracks(:,1)] *param.Xtime/param.pixel_size /1000; % Convert to pixel to overlay with raw image...
    cmap = rand(length(track_nums), 3);
    figure('Color','white')
    imagesc(max(zim,[],3)) % Fluo stack maximum projection
    hold on
    axis equal
    colormap gray
    for j = 1:length(track_nums)
        plot(xy_pix(fovTracks(:,4)==track_nums(j), 1), xy_pix(fovTracks(:,4)==track_nums(j), 2), 'LineWidth',1, 'Color', cmap(j,:));
    end
    hold off
    drawnow;
    
end

fprintf('Total number of cells: %.0f\n(including %.0f where tracking failed)\n', sum(cells_numb), skipped);


T = array2table(trackData, 'VariableNames', {'X_nm', 'Y_nm', 'Frame', 'ParticleID', 'Peak_intensity', 'FOV_number', 'CellID'});
writetable(T, [folder 'TrackData.csv']);










