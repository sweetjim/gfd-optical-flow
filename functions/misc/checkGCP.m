function hasGCP = checkGCP(state)
%CHECKGCP checks the status of the parallel-pool. 
% If a parallel pool is desired, set state == true.
if nargin==0
    state = false;
end
hasGCP = ~isempty(gcp('nocreate'));
if ~hasGCP&&state
    gcp;
end
end

