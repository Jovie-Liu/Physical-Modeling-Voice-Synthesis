function Edit(hMainFigure,parameter,pos)
uicontrol('Parent',hMainFigure, ...
        'Style','edit', ...
        'Position',pos, ...
        'Value',parameter.value, ...
        'Callback',@edit_step);
     function edit_step(edit,~)
        val = get(edit,'String');
        parameter.value = str2double(val);
    end
end