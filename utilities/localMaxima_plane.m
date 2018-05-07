% -------------------------------------------------------------------------
% [Ben] 02/22/18
% Scans through the pixels in the image array `img`, and returns a
% corresponding binary array that indicates pixels that are local maxima.
% For each pixel, considers the 3x3x3 bounding box centered on the pixel.
% If this pixel has the maximum value within the bounding box, then records
% a 1 in the same position of this pixel in corresponding binary array to
% be returned.  
% -------------------------------------------------------------------------

function BW=localMaxima_plane(Img)
[m,n,l]=size(Img);

BW=zeros(size(Img));
xRadius=2;
yRadius=2;
zRadius=2;

tempImg=zeros(m+xRadius*2, n+yRadius*2, l+zRadius*2);
tempImg(xRadius+1:m+xRadius, yRadius+1:n+yRadius, 1+zRadius:zRadius+l)=Img;

for i=xRadius+1:m+xRadius
    for j=yRadius+1:n+yRadius
        for k=1+zRadius:zRadius+l
            currRegion=tempImg(i-xRadius:i+xRadius, j-yRadius:j+yRadius, k-zRadius:k+zRadius);
            [r,c,v] = ind2sub(size(currRegion),find(currRegion == max(max(max(currRegion)))));
            if length(r)==1
                BW(i+r-2*xRadius-1, j+c-2*yRadius-1,k+v-2*zRadius-1)=1;
            end
%             for idx=1:length(r)
%                 BW(i+r(idx)-2*xRadius-1, j+c(idx)-2*yRadius-1, k+v(idx)-2*zRadius-1)=1;
%             end
        end
    end
end

end

