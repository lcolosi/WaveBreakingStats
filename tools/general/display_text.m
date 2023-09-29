function display_text(message, format)

    %%%%
    % display_text(message,format)
    %
    % Function for computing the mean and standard deviation of each pixel
    % for each track during the stable periods. Additionally outliers are 
    % identified in the mean image. 
    %
    %   Parameters
    %   ----------
    %   message : String of characters. 
    %   format  : Specifies whether the message is a title, section, body.
    %             Options include 'title', 'section', or 'body'. 
    % 
    %   Returns
    %   -------
    %   Displays message in the command window between fancy boarders for
    %   title and section and plainly fore body text.
    % 
    %%%%

    %--- Title ---% 
    if strcmp(format,'title') 

        % Print title text
        fprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',...
        '                                                                           ',...    
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
        '---------------------------------------------------------------------------',...
         message,...
        '---------------------------------------------------------------------------',...
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
        '                                                                           ');
    
    %--- Section ---%
    elseif strcmp(format,'section')

        % Print section text
        fprintf('%s\n%s\n%s\n%s\n%s\n',...
        '                                                                           ',...    
        '---------------------------------------------------------------------------',...
         message,...
        '---------------------------------------------------------------------------',...
        '                                                                           ');
    elseif strcmp(format,'body')
        
        % Print body text
        disp(message);
    end
end