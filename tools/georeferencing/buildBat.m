function buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,trackNr,res,option_plot,mask)

    %%%%
    % buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,trackNr,res,option_plot,mask)
    %
    % Function for generating Trimble georeferencing batch file for a 
    % single flight track.    
    %
    %   Parameters
    %   ---------- 
    %   dirin_img  : Path to nongeoreferenced images.    
    %   dirin_prj  : Path to Trimble project directories (the directory 
    %                containing EO, Images, Output, Project, Template).
    %   proj_lat   : Project area's mean latitude in decimal degree units. 
    %   proj_lon   : Project area's mean longitude in decimal degree units.
    %   trackNr    : Scalar track number that specifies which track will be
    %                georeferenced.
    %   res        : Resolution increment of output imagery to nearest
    %                xx cm. Highest resolution increment input is 1, trimble
    %                will start to interpolate at high resolutions and
    %                product becomes non-physical. Units: centimeters
    %  option_plot : Specifies whether to plot supplemental quality control
    %                figures. Boolean value where option_plot = 1 means 
    %                figures will be plotted. 
    %  mask        : Specifies whether the trimble project will be compiled
    %                for the equalized images with vingetting removed or 
    %                the processed image used for identifying regions of 
    %                high glint. Options of mask include: 
    %                   (1) mask = 0 : Input image is the equalized image.
    %                   (2) mask = 1 : Input image is the processed image
    %                                  for idenifying high glint.
    %                This variable will be used to specify the directories
    %                where the prj, xml, and orthomaster files will be
    %                saved. 
    % 
    %   Returns
    %   -------
    %   Generates three files:
    %
    %       (1) Applications master PRJ file: This file provides the camera
    %           specifications, trajectory information from the EO file,
    %           and specifies the datnum used (e.g., EPSG 32610). It also
    %           defines each image (i.e., path, filename, 
    %           boresight and EO attitude rotation matrix, etc.). 
    %
    %       (2) Orthomaster XML file: This file tells the project the
    %           how to export the images (i.e., where to export them, what 
    %           resoultion to export them at, what surface the image is
    %           projected onto (here it is a flat surface defined by the
    %           mea sea level).
    %
    %       (3) Batch file : This file executes the trimble project in a
    %           low level mode (specified by the flag -batch). 
    %
    %   Along with these files are a few supplementary files generated:
    %
    %       (1) Project Summary text file: For each project, a summary file
    %           is generated to output basic information about the 
    %           images processed (number of images, resolution, etc.). This
    %           file is located in the Project subdirectory. 
    %
    %       (2) Project log text file: While the images are being
    %           georeferenced a log file (with a file name similar to 
    %           om-20230719-163153_10368.txt) will be generated and updated
    %           continuously. This will provide basic information about the
    %           status of the trimble project and will tell you when the
    %           project has completed. 
    % 
    % 
    %   Instructions from Nick 
    %   ----------------------
    %   Prior to running this script, find the Trimble XML and PRJ files 
    %   from the boresight (usually on airseaserver19). From these files 
    %   you will need to make a template versions with flight-specific
    %   information not included but boresight information included 
    %   (Nick S. may have already done this, so check template directory
    %   first). You can check template versions for earlier sets of
    %   projects to get an idea of what lines you will need.
    %   
    %   Also you will need to setup all the correct directory structure
    %   for the script and batch processing to run properly. For each 
    %   flight date you will need the following folders: 
    %       (1) EO (with *EO.txt file) 
    %       (2) Images (with cropped images) 
    %       (3) Output
    %       (4) Project 
    %       (5) Templates (with *Template.xml and *Template.prj files)
    %   Note: These subdirectories are located inside the dirin_prj
    %   directory.
    % 
    %   For running the batch file, right click the file and run as
    %   Administrator.
    % 
    %   NODES explanation:
    %       (1) NODE1: Set inputs (automated here), extracts EO file,
    %                  computes average altitude and resolution, and
    %                  creates project summary.  
    %       (2) NODE2: Generates Applications Master PRJ File
    %       (3) NODE3: Generates Orthomaster XML File
    %       (4) NODE4: Generates Batch File     
    % 
    %   Additional Notes
    %   ----------------
    %   (1) The vl file refers to the raw file output from the video
    %       camera. Each flight track (also called strip) has its only 
    %       vl file and ID number. Because this code only considers a 
    %       single flight track, there will only be one vl file refered
    %       here. 
    %   (2) For specifications of the Flare camera, see the Flare 12MP 
    %       data sheet on the IO industries website: 
    %       https://www.ioindustries.com/flare-cl-downloads  
    %   (3) For more information about camera calibration, see the 
    %       matlab documentation here: 
    %       https://www.mathworks.com/help/vision/ug/camera-calibration.html
    %    
    %%%%
    
    if mask == 0
        display_text(['Flight track ', num2str(trackNr),  ': Building batch file for equalize image...'],'body')
    else 
        display_text(['Flight track ', num2str(trackNr),  ': Building batch file for high sun-glint mask...'],'body')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1st NODE: Inputs, EO file extraction, average altitude and resolution, and project summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %--- User Inputs ---% 

    % Miscellaneous constants
    msl   = egm96geoid(proj_lat,proj_lon);                                  % Mean sea level with respect to the ellipsoid at the prescribed location (units: meters; overwrite orthometric height for terrestrial values that are above the mean sea level)                
    leap = 18;                                                              % Current number of leap seconds
    % USER_ID = 'LVC';                                                      % Initials of person running trimble project (needs to be implemented) 

    % Flare video camera parameters 
    fl    = 14.149173;                                                      % Focal length of flare video camera lens (units: millimeters)
    dw_x  = 22.53;                                                          % Sensor size in the cross-track direction (units: millimeters; the width of the sensor in the x-direction (cross-track) where the pixel size of the sensor is H x V = 4096 x 3072 and the diagonal distance across the sensor is 28.1 mm; calculated from these parameter using trigonometric relations: dw_x = 28.1*cos(atan(3072/4096)))
    pix_x = 4096;                                                           % Number of pixels in the cross-track direction 

    % Directories                                  
    file_eo  = dir([dirin_prj 'EO\' '*' '_' num2str(trackNr) '_EO.txt']);   % Directory of EO file exported from Inertial Explorer for track number trackNr
    temp_xml = dir([dirin_prj 'Template\*Template.xml']);                   % Orthomaster XML template file for this particular installation / set of flights
    temp_prj = dir([dirin_prj 'Template\*Template.prj']);                   % Project template file for this particular installation /set of flights. Derived from manual boresight project file
    
    % Set Project and XML file location
    if mask == 0
        dirV=[dirin_prj 'Project\Track_' num2str(trackNr) '\'];
    else
        dirV=[dirin_prj 'Project\TrackMask_' num2str(trackNr) '\'];
    end
    if ~exist(dirV, 'dir'), mkdir(dirV); end
    dirout_prj_xml = dirV; 
    
    % Set Orthomaster output directory (where we will save georeferenced
    % images)
    if mask == 0
        dirV=[dirin_prj 'Output\Track_' num2str(trackNr) '\'];
    else 
        dirV=[dirin_prj 'Output\TrackMask_' num2str(trackNr) '\'];
    end
    if ~exist(dirV, 'dir'), mkdir(dirV); end
    dirout_ortho = dirV; 
    
    % End User Inputs. Now Retrieve image and EO file Info
    filelist = dir([dirin_img '*.tif']);
    ind = strfind(file_eo.name,'_');
    proj_str = file_eo.name(1:ind(2));                                      %'yyyymmdd_proj_' used for file naming where proj is the project acronym that was used (e.g., SMODE or TFO)
    
    % Read EO file and extract data, skipping over variable header line length
    fid = fopen([dirin_prj 'EO\' file_eo.name]);line = '';num=0;
    while ~contains(line,'(sec)')
        num=num+1;
        line = fgetl(fid);
    end
    EO_file = textscan(fid,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);

    % Move data to EO structure and clear EO_file
    EO.file = EO_file{:,2}; EO.time = EO_file{:,1};                               % File name (MASS_VIDEO_trackNumber_frameNum) and time
    EO.east = EO_file{:,3}; EO.north = EO_file{:,4}; EO.alt = EO_file{:,5};       % UTM coordinates and altitude above geoid ellipse 
    EO.omega = EO_file{:,6}; EO.phi = EO_file{:,7}; EO.kappa = EO_file{:,8};      % Angles specifying of camera orientation relative to a flat plane
    EO.roll = EO_file{:,9}; EO.pitch = EO_file{:,10}; EO.heading = EO_file{:,11}; % Attidue angles from the IMU 
    EO.lat = EO_file{:,12}; EO.lon = EO_file{:,13}; EO.alt_msl = EO_file{:,14};   % Latitude, longitude, and altitude above mean sea level
    clear EO_file

    % Extract flight line and frame number from specific filename string
    for jj=1:length(EO.time)
            file_temp = sscanf(EO.file{jj,:},'MASS_VIDEO_%d_%d'); 
            EO.line(jj,1) = sortrows(file_temp(1)); 
            EO.frame(jj,1) = file_temp(2);
    end

    % Set up a warning message if multiple flight tracks are present
    if length(unique(EO.line)) > 1
        error('Multiple flight tracks are present in the EO file! Only one track can be processed at a time.')
    end
    
    % Compute statistics for the entire track (number of frames, average altitude above MSL in meters)
    EO.frames = length(EO.file);
    EO.alt_avg = mean(EO.alt_msl);
    
    % Compute average cross track resolution
    fov_x = 2*atan2d(dw_x,2*fl);                                            % Field of view (angular width based on the equation: FOV = 2 x arctan(sensor_size/(2f)) where f is the focal length and sensor size is the pixel resolution; units: degrees)
    swath_x = 2*tand(fov_x/2)*EO.alt_avg;                                   % Flat distance viewed by the camera in the cross track direction (units: meters)
    EO.res_x = swath_x/pix_x;                                               % Average cross track resolution for the flight track (units: meters) 
    
    % Compute the nearest xx cm export size based on specified resolution
    EO.res_cm = ceil((EO.res_x * 100) / res) * res; clear swath_x 

    % Display quanlity control figure
    if option_plot == 1

        % Plotting parameters
        fontsize = 14; 
        font = 'times';

        % Generate figure
        figure('Name', ['Flight altitude and leg resolution  - Flight Track' num2str(trackNr)]);
        set(gcf,'color',[1,1,1])
    
        % Plot altitude above mean sea level and flight resolution
        plot(EO.time-EO.time(1),EO.alt_msl-msl,'.k','MarkerSize',5);
        
        % set figure attributes
        set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
        xlabel('GPS Time elapsed (sec)','fontname',font,'FontSize',fontsize)
        ylabel('Individual Photo Altitude (m)','fontname',font,'FontSize',fontsize)
        title(['Flight Track ' num2str(trackNr)],'fontname',font,'FontSize',fontsize)
        annotation('textbox',[0.13 0.62 0.3 0.3], ...
            'String',['Avg Flight Leg Res: ' num2str(EO.res_cm*10) ' mm'], ...
            'FontSize',fontsize,'FitBoxToText','on','EdgeColor','none');
        box on

    end
    
    % Write the project summary into a formatted text file
    fid = fopen([dirout_prj_xml proj_str(1:end-1) '_Summary.txt'],'w');
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf(fid,'| Flight Line |  UTC Start Time  |                    Total Files per Line: Start - End                   | Avg Alt (/100ft) | Avg Res (/5cm) | Avg Lat (deg) | Avg Lon (deg) | Avg MSL (m) |\n');
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    
    % Write vl info in the summary file
    time = datenum(proj_str(1:8),'yyyymmdd')+rem(EO.time(1),24*60*60)/(24*60*60)-leap/(24*60*60); %#ok
    line_start = datestr(time,'yyyy/mm/dd HH:MM');                                                %#ok
    fprintf(fid,'|     %2d      | %s | %5d: %30s - %30s |      %5.0fft     |    %5dcm     |    %3.4f%c   |   %4.4f%c  |    %3.1fm   |\n',unique(EO.line), line_start, EO.frames, EO.file{1,:}, EO.file{end,:}, round(EO.alt_avg*.0328084)*100, EO.res_cm, proj_lon, char(176), proj_lat, char(176), msl);
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n');
    fclose(fid);
    
    display_text('    Node 1: EO file extraction and flight summary complete.','body')
                        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2nd NODE: Applications Master PRJ File Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Store the template PRJ file in memory and find all the lines that need to be edited for automation
    fid = fopen([dirin_prj 'Template\' temp_prj.name],'r');
    ind1 = 0; ind2 = 0; ind3 = 0;
    
    while ~feof(fid)

        % Increment line counter
        ind1 = ind1 + 1; 

        % Read ith line in prj file
        template_prj{ind1} = fgetl(fid);                                    %#ok

        %--- Indices for modifying/adding new text ---%
        if contains(template_prj{ind1},'$PROJECT_NAME') || ... 
           contains(template_prj{ind1},'$STARTING_DATE') || ... 
           contains(template_prj{ind1},'$LAST_CHANGE') || ... 
           contains(template_prj{ind1},'$STRIPS') || ... 
           contains(template_prj{ind1},'$END_COORDINATE_SYSTEM') ||...
             strcmp(template_prj{ind1},'$NAVIGATION')

            % Increment line addition counter and save line number
            ind2 = ind2 + 1; 
            index_ad(ind2) = ind1;                                          %#ok
        end 

        %--- Indices for copying information from the template file ---%
        if contains(template_prj{ind1},'$ID') || ... 
           contains(template_prj{ind1},'$FOCAL_LENGTH') || ...
           contains(template_prj{ind1},'$INS_BORESIGHT_ALIGNMENT')

            % Increment line addition counter and save line number
            ind3 = ind3 + 1; 
            index_cp(ind3) = ind1;                                          %#ok
        end
    end
    fclose(fid);
    
    % Retrieve data used for some of the subsequent lines
    col_ind  = strfind(template_prj{index_cp(1)},':');                      % Find colon before Datnum and UTM zone
    crs      = template_prj{index_cp(1)}(col_ind+2:end);                    % Retrieve Datnum and UTM zone
    col_ind  = strfind(template_prj{index_cp(2)},':');                      % Find colon before project name
    cam_id   = template_prj{index_cp(2)}(col_ind+2:end);                    % Retrieve project name
    col_ind  = strfind(template_prj{index_cp(4)},':');                      % Find colon before camera focal length
    cam_fl   = template_prj{index_cp(4)}(col_ind+2:end);                    % Retrieve camera focal length
    col_ind  = strfind(template_prj{index_cp(5)},':');                      % Find colon before boresight angles
    bore_ang = split(template_prj{index_cp(5)}(col_ind+2:end));             % Retrieve camera boresight angles
    col_ind  = strfind(template_prj{index_cp(6)},':');                      % Find colon before navigation ID
    nav_id   = template_prj{index_cp(6)}(col_ind+2:end);                    % Retrieve navigation ID
    
    % Retrieve current time stamp  
    timestamp1 = datestr(now,'ddd mmm dd HH:MM:SS yyyy');                   %#ok Current timestamp to add to file
    timestamp2 = datestr(now,'HH:MM:SS dd/mm/yyyy');                        %#ok Current timestamp to add to file (slightly different format)
    
    % Replace or add lines to PRJ file as needed
    template_prj{index_ad(1)} = [template_prj{index_ad(1)} proj_str(1:end-1)]; %$PROJECT_NAME
    template_prj{index_ad(2)} = [template_prj{index_ad(2)} timestamp1];        %$STARTING_DATE
    template_prj{index_ad(3)} = [template_prj{index_ad(3)} timestamp1];        %$LAST_CHANGE
    
    % Boresight Exterior Orientation Matrix
    ai = str2double(bore_ang{end-2}); aj = str2double(bore_ang{end-1}); ak = str2double(bore_ang{end});
    bo_matrix = [cosd(aj)*cosd(ak),                            -cosd(aj)*sind(ak),                             sind(aj);...
                 cosd(ai)*sind(ak)+sind(ai)*sind(aj)*cosd(ak),  cosd(ai)*cosd(ak)-sind(ai)*sind(aj)*sind(ak), -sind(ai)*cosd(aj);...
                 sind(ai)*sind(ak)-sind(aj)*cosd(ai)*cosd(ak),  sind(ai)*cosd(ak)+cosd(ai)*sind(aj)*sind(ak),  cosd(ai)*cosd(aj)];

    %--- Add $PHOTO PRJ line information ---%
    photo(1:EO.frames,1:13) = strings;

    % Loop through frames
    for mm=1:EO.frames

        % Enter basic information
        photo{mm,1,:} = '$PHOTO'; 
        photo{mm,2,:} = sprintf('  $PHOTO_NUM : %03d%05d',EO.line(mm),EO.frame(mm));
        photo{mm,3,:} = sprintf('  $PHOTO_FILE : %s',[dirin_img filelist(mm).name]);
        photo{mm,4,:} = sprintf('  $CAMERA_ID : %s',cam_id);
        photo{mm,5,:} = sprintf('  $TERRAIN_HEIGHT : %f',msl);
        photo{mm,6,:} =         '  $IO_STAT : manual';
        photo{mm,7,:} =         '  $ORI_STAT : Initial';
        photo{mm,8,:} = sprintf('  $EXT_ORI : %s',timestamp2);
        photo{mm,9,:} = sprintf('      %10.5f   %13.5f   %13.5f   %13.5f',str2double(cam_fl),EO.east(mm),EO.north(mm),EO.alt(mm)); 

        % Calculate exterior orientation from boresight angles and EO file rounded values
        ei = round(EO.omega(mm),5); ej = round(EO.phi(mm),5); ek = round(EO.kappa(mm),5);   
        eo_matrix = [cosd(ej)*cosd(ek),                            -cosd(ej)*sind(ek),                           sind(ej);...
                     cosd(ei)*sind(ek)+sind(ei)*sind(ej)*cosd(ek),  cosd(ei)*cosd(ek)-sind(ei)*sin(ej)*sin(ek), -sind(ei)*cosd(ej);...
                     sind(ei)*sind(ek)-sind(ej)*cosd(ei)*cosd(ek),  sind(ei)*cosd(ek)+cosd(ei)*sin(ej)*sin(ek),  cosd(ei)*cosd(ej)];   
        eobo_matrix = (eo_matrix*bo_matrix)';  

        % Write to cell array
        photo{mm,10,:} = sprintf('      %15.12f      %15.12f      %15.12f',eobo_matrix(1,1),eobo_matrix(1,2),eobo_matrix(1,3));
        photo{mm,11,:} = sprintf('      %15.12f      %15.12f      %15.12f',eobo_matrix(2,1),eobo_matrix(2,2),eobo_matrix(2,3));
        photo{mm,12,:} = sprintf('      %15.12f      %15.12f      %15.12f',eobo_matrix(3,1),eobo_matrix(3,2),eobo_matrix(3,3));
        photo{mm,13,:} = '$END';  
    end 
    
    %--- Add $NAVIGATION PRJ line information ---%
    nav(1:EO.frames,1:5) = strings;

    % Loop through frames
    for qq=1:EO.frames 

        % Enter basic information
        nav{qq,1,:} = sprintf('  $PHOTO_NUM :  %03d%05d',EO.line(qq),EO.frame(qq));
        nav{qq,2,:} = sprintf('  $NAVIGATION_PAR : %s',nav_id);
        nav{qq,3,:} = sprintf('  $GPS :     %11.3f     %11.3f     %10.3f',EO.east(qq),EO.north(qq),EO.alt(qq));
        nav{qq,4,:} = sprintf('  $INS :     %10.5f     %10.5f     %10.5f',EO.omega(qq),EO.phi(qq),EO.kappa(qq));
        nav{qq,5,:} = sprintf('  $COORDINATE_SYSTEM_ID : %s',crs);      
    end 
    
    % Write new file using the template with added strip, photo, and nav data for each output resolution
    fid = fopen([dirout_prj_xml proj_str num2str(EO.res_cm) 'cm.prj'],'w');

    %--- Lines 1-25 ---%
    for aa=1:index_ad(4)
        fprintf(fid,'%s\n',template_prj{aa});
    end

    %--- Lines 26-48 ---%
    for aa=index_ad(4)+1:index_ad(5)
        fprintf(fid,'%s\n',template_prj{aa});
    end    

    %--- Write photo parameters ---%

    % Generate waitbar
    f = waitbar(0,'Please wait...');

    % Loop though photos
    for tt=1:length(photo)
    
        % Update waitbar
        if rem(tt,round(length(photo)/100,-1))==0
            waitbar(tt/length(photo),f,...
                {'Compiling photo parameters for Trimble project files'; [num2str(round((tt/length(photo))*100)) '$\%$ complete...']})
        end 
        
        % Loop through photo parameters
        for aa=1:size(photo,2)
            fprintf(fid,'%s\n',photo{tt,aa,:});
        end
    end

    % Close waitbar
    close(f)

    %--- Lines 49-93 ---%
    for aa=index_ad(5)+1:index_ad(6)
        fprintf(fid,'%s\n',template_prj{aa});
    end

    %--- Write navigation parameters ---%

    % Generate waitbar
    f = waitbar(0,'Please wait...');

    % Loop though photos
    for tt=1:length(photo)

        % Update waitbar
        if rem(tt,round(length(photo)/100,-1))==0
            waitbar(tt/length(photo),f,...
                {'Compiling navigation parameters for Trimble project files'; [num2str(round((tt/length(photo))*100)) '$\%$ complete...']})
        end 
        
        % Loop through photo parameters
        for aa=1:size(nav,2)
            fprintf(fid,'%s\n',nav{tt,aa,:});
        end
    end

    % Close waitbar
    close(f)

    %--- Line 94 ---%
    for aa=index_ad(6)+1:length(template_prj)
        fprintf(fid,'%s\n',template_prj{aa});
    end
    fclose(fid);
    
    display_text('    Node 2: Applications master PRJ file generation complete.','body')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3rd NODE: Orthomaster XML File Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Store the template XML file in memory and find all the lines that need to be edited for automation
    fid = fopen([dirin_prj 'Template\' temp_xml.name],'r');
    ind1 = 0; ind2 = 0;

    while ~feof(fid)

        % Increment line counter
        ind1 = ind1 + 1; 

        % Read ith line in prj file
        template_xml{ind1} = fgetl(fid);                                    %#ok

        %--- Indices for modifying/adding new text ---%
        if contains(template_xml{ind1},'dehPlane0') ||  ... 
           contains(template_xml{ind1},'resolution') || ... 
           contains(template_xml{ind1},'outdirPath') || ... 
           contains(template_xml{ind1},'nameMask')

            % Increment line addition counter and save line number
            ind2 = ind2 + 1; 
            index_mod(ind2) = ind1;                                         %#ok

        end    
    end
    fclose(fid);
    
    % Retrieve data used for some of the subsequent lines
    dem_ind = strfind(template_xml{index_mod(1)},'"');                      % Find first quote to insert dem value (the reference value of the flat surface you project the image onto (mean sea level))
    template_xml{index_mod(1)} = sprintf('%s%.2f%s',template_xml{index_mod(1)}(1:dem_ind),msl,template_xml{index_mod(1)}(dem_ind+1:end));
    res_ind = strfind(template_xml{index_mod(2)},'"');                      % Find first quote to insert resolution value
    template_xml{index_mod(2)} = sprintf('%s%.2f%s',template_xml{index_mod(2)}(1:res_ind(1)),EO.res_cm/100,template_xml{index_mod(2)}(res_ind(2):end));
    dir_ind = strfind(template_xml{index_mod(3)},';');                      % Find first semicolon to insert output directory
    template_xml{index_mod(3)} = sprintf('%s%s%s',template_xml{index_mod(3)}(1:dir_ind),dirout_ortho(1:end-1),template_xml{index_mod(3)}(dir_ind+1:end));
    mask_ind = strfind(template_xml{index_mod(4)},'"');                     % Find first quote to insert output name mask
    template_xml{index_mod(4)} = sprintf('%s%s%dcm_&lt;PHOTO>%s',template_xml{index_mod(4)}(1:mask_ind(1)),proj_str,EO.res_cm,template_xml{index_mod(4)}(mask_ind(2):end));   
        
    % Write new file using the template 
    fid = fopen([dirout_prj_xml proj_str num2str(EO.res_cm) 'cm.xml'],'w');    
    for aa=1:length(template_xml)
        fprintf(fid,'%s\n',template_xml{aa});
    end
    fclose(fid);

    display_text('    Node 3: Orthomaster XML file generation complete.','body')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 4th NODE: Batch File Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Write batch file to run all of the PRJ/XML files for an individual flight track 
    trimble_dir = 'C:\Program Files\Trimble\Trimble Photogrammetry 11.0\bin'; 
    prj_files = dir([dirout_prj_xml '*cm.prj']);
    fid = fopen([dirout_prj_xml proj_str(1:end-1) '.bat'],'w');
    fprintf(fid,'cd /d %s\n',trimble_dir);
    
    % Write the total number of images and the resolution to echo during georectification
    fprintf(fid,'echo %d Total Images; ',length(EO.frame));
    fprintf(fid,' at %dcm resolution; ',EO.res_cm);
    fprintf(fid,'\n');
    
    % Write the Trimble batch commands 
    fprintf(fid,'orthomaster.exe -batch -prj %s%s\n',dirout_prj_xml,prj_files.name);
    fclose(fid);

    display_text('    Node 4: Batch file generation complete.','body')

end
