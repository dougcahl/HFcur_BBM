function y = make_segments_overlap(x,N,overlap)
% function y = make_segments(x,segment_length)
%
%
%   Input:  x               = your data series
%           N               = length of segments you want
%           overlap         = percentage over lap. 50 = 50% overlap
%           keep_last_data  = 0 or 1. 0 will throw away last fraction of a
%           segment, while 1 pads it with 0s
%           if x = [1 2 3] and we want segment length of 2, overlap = 0
%           keep_last_data  = 0 will make y = [1 2]
%           keep_last_data  = 1 will make y = [1 2 ; 3 0]
%           defaults to 0 if not given
%
%   Output: y(N,n) = n segments each with a length of N
%
%
%   Extra data at the end of your data series that doesn't fit into a
%   segment will be discarded
%
%
% DLC 2014
x   = x(:);
dx  = round(N*(1-overlap/100));
m   = floor((length(x)-N)/dx) + 1;
y = complex(nan(N,m));

for i = 1:m
    y(:,i)  = x(1 + (i-1)*dx:N + (i-1)*dx);
end
