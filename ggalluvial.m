function ggalluvial(data,options)
arguments
    data
    options.cols
end



if ~isempty(options.cols)
    hexcols=rgb2hex(options.cols);
    hexcols=string(hexcols);
end