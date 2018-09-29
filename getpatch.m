function [Hp,rows,cols] = getpatch(sz,p)
% [x,y] = ind2sub(sz,p);  % 2*w+1 == the patch size
    w=5; p=p-1; y=floor(p/sz(1))+1; p=rem(p,sz(1)); x=floor(p)+1;
    rows = max(x-w,1):min(x+w,sz(1));
    cols = (max(y-w,1):min(y+w,sz(2)))';
    Hp = sub2ndx(rows,cols,sz(1));

function N = sub2ndx(rows,cols,nTotalRows)
    X = rows(ones(length(cols),1),:);
    Y = cols(:,ones(1,length(rows)));
    N = X+(Y-1)*nTotalRows;
