% -------------------------------------------------------------------------
% [Ben] 2/3/18
% Distinguishes cortical and oral apparatus BBs from BB candidates. Does
% this by looking at the skewness and kurtosis of pixel values in a window
% centered at the location of the BB candidate.
% -------------------------------------------------------------------------

% If return 1, then it's a cortical BB
% If return 0, then it's an OA BB
function result = distinguishBBClasses(img, xIdx, yIdx, patchSize)
% intensities = zeros(patchSize, patchSize);
% for i = 1:patchSize
%     for j = 1:patchSize
%         currX = xIdx+i-round(patchSize/2);
%         currY = yIdx+j-round(patchSize/2);
%         if currX > 0 && currX <= size(img,1) && currY > 0 && currY <= size(img,2)
%             intensities(i,j) = img(currX,currY);
%         end
%     end
% end
% 
% intensities_vector = reshape(intensities, [patchSize*patchSize,1]);
% intensities_vector(intensities_vector == 0) = []; % remove 0-entries

intensities = img((xIdx-patchSize/2):(xIdx+patchSize/2), (yIdx-patchSize/2):(yIdx+patchSize/2));
intensities_vector = reshape(intensities, [size(intensities, 1)*size(intensities, 2), 1]);

skew = skewness(intensities_vector);
kurt = kurtosis(intensities_vector);
% cortical BBs have a positive skew and kurtosis, while oral apparatus BBs
% have a negative skew and kurtosis due to their homogenous background. BBs
% are excluded if the sum of their skew and kurtosis is less than 1.

if (kurt + skew > 1)
    if kurt > 3
        result = 1;
    else
        result = 0;
    end
else
    result = -1;
end
end

