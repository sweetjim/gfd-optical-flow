function treecopy(tree1,tree2)
%
if nargin==0
    fig = uifigure('Name','Tree Copier Example');
    gl = uigridlayout(fig,'ColumnWidth',{'1x','1x'},'RowHeight',{'1x'});
    tree1 = uitree(gl);
    for i=1:10
        nodes = uitreenode(tree1,'Text',sprintf('Tmp%i',i));
        for j=1:ceil(rand*5)
            child1 = uitreenode(nodes,'Text',sprintf('Tmp%iChild%i',i,j));
            for k=1:ceil(rand*10)
                uitreenode(child1,'Text',sprintf('Tmp%iChild%iChild%i',i,j,k));
            end
        end
    end
    tree2 = uitree(gl);
end
%%
arrayfun(@(x,y) copyobj(flipud(x.Children),y),tree1,tree2)
end

