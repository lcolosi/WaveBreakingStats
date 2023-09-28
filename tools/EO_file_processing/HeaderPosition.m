function Counter=HeaderPosition(EOpath,String)
    
    %%%%
    % Counter = HeaderPosition(EOpath,String)
    %
    % Function for reading EO file and determine from where does the data
    % start. It is assumed that the EO file has the following data 
    % structure:
    %
    %       GPSTime Station                           Easting     Northing     H-Ell          Omega            Phi          Kappa
    %         (sec)                                       (m)          (m)       (m)          (deg)          (deg)          (deg)
    % 413939.399994 MASS_VIDEO_10_0000             312737.977  3664784.192   284.309   5.7230822161  -0.0289898081 335.9019834779
    %
    %
    %   Parameters
    %   ----------
    %   EOpath  : Path to EO text file.
    %   String  : Character string in the last row of the header
    %            (cooresponds to the row above the data).
    % 
    %   Returns
    %   -------
    %   Counter : Specifies the number of lines in the text file until 
    %             first line of data is reached (i.e. number of file lines
    %             in the header).           
    %
    %%%%

    % Open text file for reading access
    fid=fopen(EOpath);

    % Read first line of EO file
    tline = fgetl(fid);

    % Set counter
    Counter=1;

    % Loop through lines in text file and check if line contains '(sec)'
    % string
    while ~contains(tline,String)
        % Obtain ith line's character string
        tline = fgetl(fid);
        % Step forward counter
        Counter=Counter+1;
    end

    % Close file when first line with data is reached     
    fclose(fid);

end