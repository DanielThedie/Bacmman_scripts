
close all
clc

folder = '/media/daniel/HDD Daniel/Daniel Thédié/Tracking/210820/DT2_100pc/';
fileName = 'TrackData.csv';


%% End of input

% Load data
trackData = readtable([folder fileName]);

% Compute Mean Jump Distance
trackLen = groupcounts(trackData.ParticleID); % Track length = number of localisations in the track
track_nums = unique(trackData.ParticleID); % List of all particle IDs

mjd = zeros(length(trackLen),1);
track_centre = zeros(length(trackLen), 2);
jd_list = [];
for i = 1:length(trackLen) % For each track
    track_centre(i,:) = [mean(trackData{trackData.ParticleID==track_nums(i),1}) mean(trackData{trackData.ParticleID==track_nums(i),2})];
    
    D = pdist([trackData{trackData.ParticleID==track_nums(i),1} trackData{trackData.ParticleID==track_nums(i),2}]);
    Z = squareform(D);
    ind = 1:sum(trackData.ParticleID==track_nums(i));
    dt = triu(ind-ind'); % Make a table of distance between different items
    jd = Z(dt == 1);
    
    jd_list = [jd_list; jd];
    
    mjd(i) = mean(jd);
end
mjd_n = trackLen -1; % Number of jumps per track = number of localisations -1

% Weight MJD by number of jumps per track
h_data = zeros(sum(mjd_n), 1);
c = 1;
for j = 1:length(mjd)
    h_data(c:c+mjd_n(j)) = mjd(j);
    c = c + mjd_n(j);
end

% Make histogram
figure('Color', 'white')
histogram(h_data, 'Normalization', 'probability')
xlabel('MJD (nm)')
ylabel('PDF')
















