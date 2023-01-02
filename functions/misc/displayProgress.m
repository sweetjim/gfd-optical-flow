function varargout = displayProgress(string,val,start_val,end_val,varargin)

if nargout==0
    if val==start_val
        fprintf('\r%s:\t',string)
    end
end
abdiff = @(i) diff([i end_val]);
frac = round((abdiff(start_val)-abdiff(val))/abdiff(start_val)*100);

if nargout>0
    varargout{1} = sprintf('%s = %i%%',string,frac);
    return
end

if frac == 0; frac=1;end
if val>start_val
    for j=0:log10(frac*1e3)
        fprintf('\b'); 
    end
end
lines = fprintf(' %i %%', frac);
pause(.001); 

if val==end_val
    if nargin>4
        if strcmp(varargin{1},'delete')
            fprintf(repmat('\b',1, lines+length(string)+1))
        else
            fprintf('\n');
        end
    else
        fprintf('\n');
    end
end
end

