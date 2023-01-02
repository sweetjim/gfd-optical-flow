function hFound = pointer2object(figurehandle,objecttags)
%POINTER2OBJECT Summary of this function goes here
%   Detailed explanation goes here
pos         = figurehandle.CurrentPoint;
hFound      = [];
tmp         = findobj(figurehandle);

for n=1:length(objecttags)
    for j=1:length(tmp)
        try %#ok<TRYNC>
            if strcmp(tmp(j).Tag,objecttags{n})
                idx = j;
                break
            end
        end
    end
    
    hTag = tmp(idx);
    hPos = getpixelposition(hTag,true);  % recursive to find absoulute position
    if pos(1)>hPos(1) && pos(1)<hPos(1)+hPos(3) && pos(2)>hPos(2) && pos(2)<hPos(2)+hPos(4)
        hFound = hTag;
        break
    end
end
% if isempty(hFound)
%     return
% end

end

