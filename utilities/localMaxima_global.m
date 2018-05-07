% -------------------------------------------------------------------------
% [Ben] 12/06/17
% Basically a wrapper for imregionalmax. Identifies regional maxima - 
% connected components of pixels with constant intensity value, whose
% boundary pixels all have lower value. Returns image of same dimensions
% except that pixels part regional maxima have their value set to 1, all
% else is 0. 
% -------------------------------------------------------------------------

function BW = localMaxima_global(Img)
BW = imregionalmax(Img);
end
