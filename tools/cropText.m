function An=cropText(A)

    %%%%
    % An = cropText(A)
    %
    % Function for cutting out parts of the header that were not removed 
    % when importing data (look for '(sec)' in the string in the first 
    % column of A.textdata). The data field in the A structure is not 
    % contaminated with header strings (A.data contains the position and 
    % attitude data). Only the textdata is cropped (A.textdata contains the
    % GPSTime station and contains header strings). 
    %
    %   Parameters
    %   ----------
    %   A : Aircraft trajectory and atitude data with uncropped GPSTime 
    %       station cell array (A.textdata).
    % 
    %   Returns
    %   -------
    %   A : Cropped GPSTime station cell array.            
    %
    %%%%

    % Set counter
    Counter=2;
    
    % Obtain first line in the cell array
    tline=A.textdata(1,1);
    
    % Loop through rows in cell array until row contains '(sec)' characters 
    % (signifying the end of the header)
    while ~contains(tline,'(sec)')

        % Set next row character string
        tline=A.textdata(Counter,1);
        
        % Step forward counter
        Counter=Counter+1;
    end
    
    % Crop GPSTime Station data
    An=A.textdata(Counter:end,:);