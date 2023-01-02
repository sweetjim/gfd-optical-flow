function out = within(x,xrange)
%% Within
% Returns a logical argument that depends on whether 'x' falls within the
% range of 'xrange'.
% 
% Example:
% >> within(0,[-1 1])
% >> ans =
% >>    logical
% >>        1
% --------------------
% >> within(5.1,0:5)
% >> ans =
% >>    logical
% >>        0
%%
x1 = xrange(1);
x2 = xrange(end);
out = (x>=x1&x<=x2);
if diff(xrange)<0
out = (x<=x1&x>=x2);
end
end

