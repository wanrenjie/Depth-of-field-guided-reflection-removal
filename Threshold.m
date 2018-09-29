function Th = Threshold(img)
    [height,width] = size(img);
    Th = (1/(height*width))*sum(img(:));
end
