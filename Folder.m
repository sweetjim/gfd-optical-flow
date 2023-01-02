classdef Folder
    %% FOLDER is a class that constructs a multi-level variable containing all nested subfolders of the defined root folder
    %%
    properties
        Path        % Path relative to the input
        Level       % Level of the folder relative to the input
        Label       % Label of the current folder
        Parent      % Parent of the current folder [Folder class]
        Children    % Children of the current folder [cell array of Folder classes]
        Files       % Files of the current folder [DIR structure array]
    end    
    properties (Hidden)
       allfolders               % List of all folders nested under input folders
       isTree       = false     % Boolean for tree object
       useFileExt   = {}        % File extensions go here
       useChecker   = false     % Boolean for check-tree object
    end
    
    methods
        function obj = Folder(path,varargin)
            %% CONSTRUCTOR
            %% BUILD DIRECTORY
            parent      = {};
            % persistent iterations % DEBUG TOOL
            parseInput(varargin)
            
            % Remove '..' from paths
            if contains(path,'\..')
                [~,tmp] = fileattrib(path);
                path =  tmp.Name;
            end

            obj.Path    = path;
            list        = split(genpath(path),';');
            list        = list(1:end-1);
            
            % ITERATIONS (DEBUGGING)
            if isempty(parent)
                obj.Level   = 1;
            else
                obj.Level   = count(path,'\')+1;
            end
            
            % LABEL
            obj.Label   = split(path,'\');
            obj.Label   = obj.Label{obj.Level};
            
            % PARENT & FILES
            if obj.Level==1
                parent = obj;
                tmp = dir(obj.Path);
                obj.allfolders = list;
            else
                obj.Parent = parent;
                tmp = dir(fullfile(parent.Path,obj.Label));
            end
            tmp = tmp(3:end);
            
            if isempty(tmp)
                obj.Files       = {};
            else
                obj.Files       = tmp(~[tmp.isdir]);
            end
            
            % CHILDREN
            try
                if length(list)==1
                    []; %#ok<VUNUS> % DEBUG POINT
                elseif length(list)>1
                    for i=1:count(list(2),'\')
                        list = extractAfter(list,'\');
                    end
                    tmpChildren     = list(2:end);
                    omit    = logical(zeros(1,length(tmpChildren))); %#ok<LOGL>
                    tmp     = extractBefore(tmpChildren,'\');
                    for i=1:length(tmp)
                        omit(i)=isempty(tmp{i});
                    end
                    tmpChildren = tmpChildren(omit);
                end
                % RECURSION
                for i=1:length(tmpChildren)
                    tmp             = Folder(...
                        fullfile(obj.Path,tmpChildren{i}),...
                        'parent',obj);
                    
                    if i==1 % root
                        if length(tmpChildren)==1
                            folderChildren = tmp;                            
                        else % root with children
                            folderChildren = {tmp};
                        end
                    else    % combine children folders
                        folderChildren  = vertcat(folderChildren,{tmp}); %#ok<AGROW>
                    end
                    
                end
                obj.Children = folderChildren;
            catch
                obj.Children    = {};
            end
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'parent'
                            parent = namevalue;
                    end
                    m = m+1;
                    if m>length(items);break;end
                end
                function out = namevalue
                    out = items{m+1};
                    m   = m+1;
                end
            end
        end
        function maketree(obj,varargin)
            %% BUILD A UITREE FROM ROOT DIRECTORY
            % 
            %% Preamble
            handle      = [];
            makefiles   = false;
            allowext    = {};
            icons       = true;
            update      = false;
            rename      = false;
            renamestruct=struct('RootNode',[],'NodeData',[]);
            childpath   = {};
            replacepath = {};
            checkbox    = false;
            checkanalysis = false;
            parseInput(varargin)
            
            obj.isTree      = true;
            obj.useFileExt  = allowext;
            obj.useChecker  = checkanalysis;
            %
            if icons
                iconpath = fullfile('icons/folder_win10.png');
                fileicon = fullfile('icons/icon-19.png');
            else
                iconpath = '';
                fileicon = '';
            end
            
            if isempty(handle)
                tree = uitree;
                handle = tree;
            end
            
            if update
                if ~isempty(replacepath)
                    replace = childpath{ismember(childpath,setxor(childpath,replacepath))};
                    with    = [];
                end
                
                folderhandle = findfolder(obj,handle,childpath);
                % to add..
                return
            elseif rename
                %% Rename folder protocol
                rootnode    = renamestruct.RootNode;
                omit        = findall(rootnode);
                delete(omit(omit~=rootnode))
                rootnode.NodeData = renamestruct.NodeData;
                makeChildNodes(rootnode.NodeData,rootnode)
                sortNodes(rootnode)
                %%
                return
            end
            %% GENERATION
            delete(handle.Children)
            parent  = handle;
            if parent.Editable
                parent.NodeTextChangedFcn = @treeRenameNode;
            end
            makeChildNodes(obj,parent)
            sortNodes(handle)
            
            %% Functions
            function sortNodes(handle)
                % Sort the tree by numeric/alphabetic order
                if isempty(handle.Children)
                    return
                end
                children        = {handle.Children.Text}';
                if numel(children)==1
                    return
                end
                [~,sortedidx]   = sort(cellfun(@lower,children,'UniformOutput',false));
                unsortedidx     = (1:max(sortedidx))';
                mapper          = [unsortedidx sortedidx];
                for ii=1:length(mapper)
                    if mapper(ii,1)~=mapper(ii,2)
                        move(handle.Children(mapper(ii,2)),handle.Children(mapper(ii,1)))
                        children        = {handle.Children.Text}';
                        [~,sortedidx]   = sort(cellfun(@lower,children,'UniformOutput',false));
                        unsortedidx     = (1:max(sortedidx))';
                        mapper          = [unsortedidx sortedidx];
                    end
                end
            end
            %% NODE: CREATION / RECURSION
            function makeChildNodes(node,parent)
                arrayfun(...
                    @(x) createNode(parent,x),...
                    node)
                %%
                function createNode(parent,child)
                    hasChildren = ~isempty(child.Children);
                    hasFiles    = ~isempty(child.Files);
                    if hasChildren
                        for i=1:numel(child.Children)
                            newnode = uitreenode(parent,...
                                'Text',child.Children{i}.Label,...
                                'NodeData',child.Children{i},...
                                'UserData',child.Children{i}.Path,...
                                'Icon',iconpath);
                            makeChildNodes(child.Children{i},newnode)
                            treenodeCM(newnode)
                        end
                    end
                    if hasFiles&&makefiles
                        data = child.Files;
                        for j=1:length(data)
                            if ~isempty(regexp(data(j).name,allowext)')
                                newnode = uitreenode(parent,...
                                    'Text',data(j).name,...
                                    'Icon',fileicon,...
                                    'UserData',child.Files(j));
                                treenodeCM(newnode)
                            end
                        end
                    end
                end

            end
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'handle'
                            handle = namevalue;
                        case 'files'
                            makefiles = true;
                        case 'filter'
                            allowext = namevalue;
                        case 'icons'
                            state = namevalue;
                            switch state
                                case 'on'
                                    icons = true;
                                case 'off'
                                    icons = false;
                            end
                        case 'rename'
                            rename      = true;
                            renamestruct.RootNode = namevalue;
                            renamestruct.NodeData = namevalue;
                        case 'update'
                            update      = true;
                            childpath   = namevalue;
                            if numel(childpath)>1
                                replacepath = childpath{2};
                                childpath   = childpath{1};
                            end
                        case {'select','check','checkbox'}
                            checkbox    = true;
                        case 'analysis'
                            checkanalysis = true;
                    end
                    m = m+1;
                    if m>length(items);break;end
                end
                function out = namevalue
                    out = items{m+1};
                    m   = m+1;
                end
            end
        end
    end
    
    methods (Hidden)
        function folderhandle = findfolder(obj,parent,folder2find)
            foldersinlevel = obj;
            search(foldersinlevel)
            
            function search(folderlevel)
                for i=1:numel(folderlevel)
                    [];
                    hasChildren = ~isempty(folderlevel.Children);
                    if hasChildren
                        [];
                        folderpaths = cellfun(@(x) x.Path,folderlevel.Children,'UniformOutput',false);
                    else
                        folderpaths = folderlevel.Path;
                    end
                    
                    if any(contains(folderpaths,folder2find))
                        [];
                    else
                        if isempty(folderlevel.Children{i})
                            continue
                        end
                        search(folderlevel.Children{i})
                    end
                end
            end
        end
    end
end

