function [spots, cells] = Get_spots(segFile_param, position)


info = h5info([segFile_param.folder segFile_param.fileName]);
dataset = info.Groups(1).Groups.Groups(position).Name;
try
    dataSpots = h5read([segFile_param.folder segFile_param.fileName],...
        [dataset '/' segFile_param.spotsFeature]);
    dataCells = h5read([segFile_param.folder segFile_param.fileName],...
        [dataset '/' segFile_param.cellsFeature]);
catch
    error(['Could not find position ' num2str(position)])
end

spots = struct;
spots.Connectivity = 8;
[sz1, sz2, ~] = size(dataSpots); % For compatibility: "size(data, [1, 2])" introduced in R2019b
spots.ImageSize = [sz1 sz2];
spots.NumObjects = max(dataSpots, [], 'all');
spots.PixelIdxList = cell(1, spots.NumObjects);

for i=1:spots.NumObjects
    spots.PixelIdxList{i} = find(dataSpots == i);
end


cells = struct;
cells.Connectivity = 8;
[sz1, sz2, ~] = size(dataCells); % For compatibility: "size(data, [1, 2])" introduced in R2019b
cells.ImageSize = [sz1 sz2];
cells.NumObjects = max(dataCells, [], 'all');
cells.PixelIdxList = cell(1, cells.NumObjects);

for i=1:cells.NumObjects
    cells.PixelIdxList{i} = find(dataCells == i);
end


end