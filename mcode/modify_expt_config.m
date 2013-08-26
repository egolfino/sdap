function modify_expt_config(exptConfigFN, item, value)
ctxt = textread(exptConfigFN, '%s', 'delimiter', '\n');

ctxt1 = {};
for i1 = 1 : numel(ctxt)
    tline = deblank(ctxt{i1});
    if isempty(tline) || isequal(tline(1), '%');
        ctxt1{end + 1} = tline;
    else
        if ~isempty(strfind(tline, '%'))
            idxpct = strfind(tline, '%');
            idxpct = idxpct(1);
            
            lpt1 = tline(1 : idxpct - 1);
            lpt2 = tline(idxpct : end);
        else
            lpt1 = tline;
            lpt2 = '';
        end
        
        titems = splitstring(lpt1);
        if length(titems) == 2 && isequal(titems{1}, item)
            lpt1 = sprintf('%s %f', item, value);
        end
        
        new_line = sprintf('%s %s', lpt1, lpt2);
        ctxt1{end + 1} = new_line;
    end
end

f = fopen(exptConfigFN, 'wt');
for i1 = 1 : numel(ctxt1)
    fprintf(f, '%s\n', ctxt1{i1});
end

return