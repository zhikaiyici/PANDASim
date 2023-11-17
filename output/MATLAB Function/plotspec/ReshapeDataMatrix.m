function spectraData = ReshapeDataMatrix(arraySize, data)

eventNum = size(data, 1) ./ arraySize;
% spectraDataCell = mat2cell(spectradata, ones(eventnum, 1) .* arraysize);
data = data';
data = reshape(data,arraySize, arraySize,eventNum);
spectraData = permute(data, [2,1,3]);
end
