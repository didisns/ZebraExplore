function y = smoothGab(x, win, varargin)
%This function smoothes a vector or a matrix by averaging in a window with
%a width of "win" points.

%INPUT
%     x = data to be smoothed. 
%         if x is a matrix, the function treats columns as vectors and smoothes them.
%             
%     win = number of points for the smoothening
%     
%     varargin = name - value pairs.
%           'side';     'center' (default) = smoothening is performed over a window centered on every point.
%                       'left' = smoothening is performed over a window such that every point is the right extreme.
%                       'right' = smoothening is performed over a window such that point is the left extreme.
%
%           'weight';    '' (default) =
%                       '' =
%OUTPUT
%     y = smoothed data.

narg = size(varargin);
if mod(narg,2)~=0
    error('uncorrect name-value pairs')
end
%controlla side; per ora center
halfWin = floor(win/2);
xDim = size(x);
if min(xDim) == 1 %if x is a row vector, make it a column
    x = x';
end
xDim = size(x);
y = zeros(xDim);
for j = 1:xDim(1)
    
    left = j-halfWin;
    if left <1
        left = 1;
    end

    right = j+halfWin;
    if right > xDim(1)
        right = xDim(1);
    end

    y(j,:) = mean(x(left:right,:),1);
end
return
            
            
        