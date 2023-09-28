function utc_time = gpssw2utcdn(gps_time,date)

    %%%%
    % utc_time = gpssw2utcdn(gps_time,date)
    %
    % Function for converting GPS time in seconds of the week to UTC datenum.
    % Here, weeks since 0 hours (midnight) Sunday 6 January 1980 are not
    % included as an input so the starting date in UTC must be provided as 
    % an argument. Note that seconds of the week reset to zero on 0 hours
    % (midnight) Sunday. 
    %
    %   Parameters
    %   ----------
    %   gps_time : GPS time array in seconds of the week. 
    %   date     : UTC starting date for gps time. 
    %
    %   Returns
    %   -------
    %   utc_time : UTC datenum array.  
    %
    %   Issues
    %   ------
    %   GPS time will reset to zero at the end of the week (0 hours
    %   (midnight) Sunday). This needs to be accounted for when computing
    %   seconds of the day. I need to find a way to know when gps_time
    %   jumps to zero and to add the total seconds in a day (86400) to all
    %   the time steps after the zero reset. 
    %
    %%%%

    % Set variables for time conversion 
    T_day = 86400;                                                          % Number of seconds in a day
    
    % Compute the number of days since the beginning of the week for the
    % first gps time value  
    ndays = floor(gps_time(1)/T_day);

    % Remove seconds in ndays from the gps time to obtain seconds of the
    % day 
    seconds_of_the_day = gps_time - ndays*T_day; 

    % Convert gps time into UTC time (without accounting for leap seconds)
    utc_time_wo_ls = datenum(date,'yyyymmdd') + seconds_of_the_day*(1/T_day); 

    % Add leap seconds 
    utc_time = gps2utc(utc_time_wo_ls);

end





