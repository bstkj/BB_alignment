% -------------------------------------------------------------------------
% [Ben] 02/22/18
% Scans through the pixels in the image array `img`, and returns a
% corresponding binary array that indicates pixels that are local maxima.
% For each pixel, considers the 3x3x3 bounding box centered on the pixel.
% If this pixel has the maximum value within the bounding box, then records
% a 1 in the same position of this pixel in corresponding binary array to
% be returned.  
% -------------------------------------------------------------------------

function BW=localMaxima_global2(Img)
% BW = imregionalmax(Img);
[m,n,k]=size(Img);
BW=zeros(size(Img));
% (a,b,c) is the middle point of the box
for a=1:m
    for b=1:n
        for c=1:k
            tempMidInt=Img(a,b,c);
            tempBox=zeros(3,3,3);
            for da=-1:1
                if a+da>=1 && a+da<=m
                    for db=-1:1
                        if b+db>=1 && b+db<=n
                            for dc=-1:1
                                if c+dc>=1 && c+dc<=k
                                    tempBox(da+2, db+2, dc+2)=Img(a+da, b+db, c+dc);
                                end
                            end
                        end
                    end
                end
            end
            tempBox=reshape(tempBox, [3*3*3, 1]);
            tempBox(tempBox==0)=[];
            if tempMidInt==max(tempBox)
                BW(a,b,c)=1;
            end
        end
    end
end

           
end