function time = set_sliderlimits(app)
info                    = ncinfo(app.openExpt);
dim                     = contains({info.Dimensions.Name},'time');
time                    = info.Dimensions(dim).Length;
end

