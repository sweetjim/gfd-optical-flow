function treenodeCM(parent)
%% Construct context menu for tree nodes
% Has the following attributes:
% 1. Show in explorer
% 2. Rename
% 3. New folder

%% Figure parent
figparent = parent;
while true
    try
       figparent = figparent.Parent;
       if contains(class(figparent),'Figure')
          break 
       end
    catch
       return
    end
end

%% Creation
NodeContextMenu = uicontextmenu(figparent);

%% Create OpeninexplorerMenu
OpeninexplorerMenu = uimenu(NodeContextMenu);
OpeninexplorerMenu.MenuSelectedFcn = @treenodeCMCallback;
OpeninexplorerMenu.Text = 'Open in explorer';
OpeninexplorerMenu.Tag = 'command-explorer';

%% Create CopyImageProtocolMenu
CopyImageProtocol = uimenu(NodeContextMenu);
CopyImageProtocol.MenuSelectedFcn = @treenodeCMCallback;
CopyImageProtocol.Text = 'Copy images';
CopyImageProtocol.Tag = 'command-copyprotocol';

%% Create ROItransferMenu
ROItransfer = uimenu(NodeContextMenu);
ROItransfer.MenuSelectedFcn = @treenodeCMCallback;
ROItransfer.Text = 'ROI transfer';
ROItransfer.Tag  = 'command-ROI';

%% Create controlsTransferMenu
controlsTransferB = uimenu(NodeContextMenu);
controlsTransferB.MenuSelectedFcn = @treenodeCMCallback;
controlsTransferB.Text = 'Controls transfer';
controlsTransferB.Tag  = 'command-controls';

%% Create RenameMenu
% RenameMenu = uimenu(NodeContextMenu);
% RenameMenu.MenuSelectedFcn = @treenodeCMCallback;
% RenameMenu.Text = 'Rename';
% RenameMenu.Tag = 'command-rename';

%% Create NodeNewfolderMenu
% NodeNewfolderMenu = uimenu(NodeContextMenu);
% NodeNewfolderMenu.MenuSelectedFcn = @treenodeCMCallback;
% NodeNewfolderMenu.Text = 'New folder';
% NodeNewfolderMenu.Tag = 'command-mkdir';

%% UNFINISHED
% 
% % Create NodeNewfolderMenu
% NodeNewfolderMenu = uimenu(NodeContextMenu);
% NodeNewfolderMenu.MenuSelectedFcn = @treenodeCMCallback;
% NodeNewfolderMenu.Text = 'New folder';
% NodeNewfolderMenu.Tag = 'command-mkdir';

%% Assign NodeContextMenu
parent.ContextMenu = NodeContextMenu;

    function treenodeCMCallback(source,event)
        tag = event.Source.Tag;
        command = extractAfter(tag,'-');
        
        switch command
            case 'explorer'
                winopen(parent.NodeData.Path)
            case 'copyprotocol'
                imageCopier(parent.NodeData)
            case 'ROI'
                roiTransfer(parent.NodeData)
            case 'controls'
                controlsTransfer(parent.NodeData)
            case 'rename'
%                 th = ancestor(parent,'uitree');
%                 parent
%                 treeRenameNode(th)
            case 'mkdir'
                
                
        end
    end
end
