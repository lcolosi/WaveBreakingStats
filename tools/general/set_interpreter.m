function set_interpreter(interpreter)
    
    %%%%
    % set_interpreter(interpreter)
    %
    % Function for setting the default interpreter for all text displayed 
    % in plots.  
    %
    %   Parameters
    %   ----------
    %   Interpreter : Specifies the text interpreter of the text. Options
    %                 include: 
    %                   'none'  : Display literal characters.
    %                   'tex'   : Interpret characters using a subset of TeX markup.
    %                   'latex' : Interpret characters using LaTeX markup.      
    % 
    %   Returns
    %   -------
    %   None
    %%%%

    % Obtain field names of the graphics root objective
    list_factory = fieldnames(get(groot,'factory'));
    
    % Find field names that associated with the interpreter
    index_interpreter = find(contains(list_factory,'Interpreter'));
    
    % Loop through interpreter field names
    for i = 1:length(index_interpreter)
        
        % Set the default field name of the ith groot interpreter property 
        default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
        
        % Set the default interpreter of the ith groot interpreter property
        set(groot, default_name, interpreter);

    end