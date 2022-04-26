function read_struct = AddDefault(read_struct, default_struct)
    fields = fieldnames(default_struct);
    for f = 1:numel(fields)
        f_name = fields{f};
        if ~isfield(read_struct, f_name)
             read_struct.(f_name) = default_struct.(f_name);
        end
    end
end