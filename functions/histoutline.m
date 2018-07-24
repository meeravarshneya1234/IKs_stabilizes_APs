function h = histoutline(xx, num, varargin)
% Plot a histogram outline.
% h = HISTOUTLINE(xx, num, varargin) Plot a histogram outline.
%
% HISTOUTLINE uses hist to do most of the work:
%
% HISTOUTLINE(...) produces a plot of the outline of the histogram of the
% results. 
%
% Example: 
%     data = randn(100, 1); 
%     histoutline(data, 50);
%
% See also HIST, HISTC, MODE.    
%
% Matt Foster <ee1mpf@bath.ac.uk>

if nargin < 1
    error('No input given');
end

% Default to standard number of histogram bins.
if nargin < 2
    num = 10;
end

[n,x] = hist(xx, num);
h = stairs(x, n, varargin{:});
