%% SEQUENCER is a script that sequentially writes 
% velocity field measurements to .nc files. 
% 
%%
% change these paths or use commented-out UIGET lines
root        = '..\..\Data\SST_FLIR_corrected_netcdf\LRA_RotatingRidge\data';
savepath    = '..\..\Data\SST_FLIR_corrected_netcdf\LRA_RotatingRidge\uv_fields';   % ensure that the results folder exists!!
% root        = uigetdir(cd,'Select the root file location');
% saveloc     = uigetdir(cd,'Select the save destination folder');
files       = dir(fullfile(root,'*.nc'));
readpath    = @(x) fullfile(x.folder,x.name);
tfarray     = zeros(1,numel(files));

%
useparallel = false;        % parallel processing option
maskannulus = true;         % output has NaN mask if true

maskstr = '';
if maskannulus
    maskstr = 'mask';
end


% loop through experiments and write .nc files
for i=1:numel(files)
    data = exptdata(readpath(files(i)));            % initialize class
    data.loadImage(1);                              % initialize sizes
    t60s = roundto(data.times,60);                  % sample every 60s
    loop = 1:2;%:numel(t60s);                           % parfor-loop sequencer
    u   = zeros([numel(loop) size(data.image)]);    % initialize u,v
    v   = u;                                        %   use (t,x,z) indexing for performance
    % U   = u; Vor = u;                             % uncomment for U,vor 
    tic                                             % start stopwatch
    
    % loop through images and calculate velocities
    if useparallel&&checkGCP(useparallel)
        parfor j=loop 
            output = data.getVelocity(loop(j),'method','slow',maskstr); %#ok<PFBNS>
            u(j,:,:) = output.u;
            v(j,:,:) = output.v;
            %U(j,:,:) = out.U;                      % uncomment for U,vor
            %Vor(j,:,:) = out.vor;
        end
    else
        for j=loop
            output = data.getVelocity(loop(j),'method','slow',maskstr);
            u(j,:,:) = output.u;
            v(j,:,:) = output.v;
            %U(j,:,:) = out.U;
            %Vor(j,:,:) = out.vor;
        end
    end

    u = permute(u,[2 3 1]);                         % re-order arrays (x,z,t)
    v = permute(v,[2 3 1]);
    %U = permute(U,[2 3 1]);
    %Vor = permute(Vor,[2 3 1]);
    
    tfarray(i) = toc;                               % save elapsed time

    % write to netcdf
    filenameext = sprintf('%s_uvfield.nc',extractBefore(data.folder.name,'.nc'));
    savepath    = fullfile(saveloc,filenameext);
    if isfile(savepath)
        warning('File already exists\n%s\n',savepath)
        delete(savepath) % ?? keep if you want to overwrite
    end
    ncfileProtocol('u',u)
    ncfileProtocol('v',v)

    displayProgress('Progress',i,1,numel(files)) % print to command line percentage complete
end
