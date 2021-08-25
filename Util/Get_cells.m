function segmentation = Get_cells(segFile_param, position)

info = h5info([segFile_param.folder segFile_param.fileName]);
dataset = info.Groups(1).Groups.Groups(position).Name;
try
    data = h5read([segFile_param.folder segFile_param.fileName],...
        [dataset '/' segFile_param.cellsFeature]);
catch
    error(['Could not find position ' num2str(position)])
end

segmentation = struct;
segmentation.Connectivity = 8;
[sz1, sz2, ~] = size(data); % For compatibility: "size(data, [1, 2])" introduced in R2019b
segmentation.ImageSize = [sz1 sz2];
segmentation.NumObjects = max(data, [], 'all');
segmentation.PixelIdxList = cell(1, segmentation.NumObjects);

for i=1:segmentation.NumObjects
    segmentation.PixelIdxList{i} = find(data == i);
end

end