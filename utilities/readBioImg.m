% -------------------------------------------------------------------------
% [Ben] (12/06/17)
% Seems to be extracting pixel data from 3D image, and then reducing the
% dimensions of the representation by sampling every other plane and then
% reconstructing.
% -------------------------------------------------------------------------

function I=readBioImg(fileName, channel)
data=bfopen(fileName);
series1 = data{1, 1};
series1_planeCount = size(series1, 1);

series1_plane1=double(series1{1, 1}); % change type: uint16 to double
[m,n]=size(series1_plane1); % get size of single plane

% intensities in data_3d are the same as original dataset
I=zeros(m,n,series1_planeCount/2); % only consider one channel

for plane=channel:2:series1_planeCount
    currPlane = double(series1{plane, 1});
    for i=1:m
        for j=1:n
            I(i,j,floor((plane+1)/2))=currPlane(i,j);
        end
    end
end
end

