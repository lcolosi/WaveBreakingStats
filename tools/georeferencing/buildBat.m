%% 1st NODE: Inputs, EO file extraction, average altitude and resolution, and project summary

function buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,trackNr,camName)


%User Input File Information for a single flight. Each flight should have subdirectories EO (with *EO.txt file), Templates (with *Template.xml and *Template.prj files), Project, and Output
% dirin_img = 'D:\MASS\RAW\LCDRI_2017\VIDEO\TIFF_16BIT_EXPORT\20170323\'; %directory of all tif imagery to use for this single flight
% dirin_prj = 'D:\MASS\Processed\LCDRI_2017\VIDEO\20170323\';

% dirin_img = 'D:\MASS\RAW\LCDRI_2017\VIDEO\TIFF_16BIT_EXPORT\20170323N\'; %directory of all tif imagery to use for this single flight
% dirin_prj = 'D:\MASS\Processed\LCDRI_2017\VIDEO\20170323N\';

%User Input Constants
% proj_lat = 32; %project area latitude decimal deg
% proj_lon = -122; %project area longitude decimal deg
msl = egm96geoid(proj_lat,proj_lon); %msl with respect to the ellipsoid at the prescribed location
%msl = 38; %overwrite orthometric height for terrestrial values that are above the mean sea level
%cam = 'jai'; %current camera options are 'jai', 'flare', 'sc6000', 'sc6700', 'd810', 'doppvis'
cam = camName;

%Directory Info
eoName=['_' num2str(trackNr) '_EO.txt'];
file_eo =        dir([dirin_prj 'Trimble_EO\' '*' '_' num2str(trackNr) '_EO.txt']); %Full directory of EO file exported from Inertial Explorer for all image frames
temp_xml =       dir([dirin_prj 'Template\*Template.xml']); %Orthomaster XML template file for this particular installation / set of flights
temp_prj =       dir([dirin_prj 'Template\*Template.prj']); %Project template file for this particular installation /set of flights. Derived from manual boresight project file

dirV=[dirin_prj 'Project\Track_' num2str(trackNr) '\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end
dirout_prj_xml =     dirV; %Project and XML file location

dirV=[dirin_prj 'Output\Track_' num2str(trackNr) '\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end
dirout_ortho =       dirV; %Orthomaster output directory

%End User Inputs. Now Retrieve File Info
filelist = dir([dirin_img '*.tif']);
ind = strfind(file_eo.name,'_');
proj_str = file_eo.name(1:ind(2)); %'yyyymmdd_proj_' used for file naming where proj is the project acronym that was used

%Read EO file and extract data, skipping over variable header line length, and sort rows into order by vl/frame
fid = fopen([dirin_prj 'Trimble_EO\' file_eo.name]);line = '';num=0;
while ~contains(line,'(sec)')
    num=num+1;
    line = fgetl(fid);
end
EO_file = textscan(fid,'%f MASS_VIDEO_%d_%d %f %f %f %f %f %f');
fclose(fid);
[EO.vl,index] = sortrows(EO_file{:,2}); EO.time = EO_file{:,1}(index); EO.frame = EO_file{:,3}(index);
EO.east = EO_file{:,4}(index); EO.north = EO_file{:,5}(index); EO.alt = EO_file{:,6}(index);
EO.omega = EO_file{:,7}(index); EO.phi = EO_file{:,8}(index); EO.kappa = EO_file{:,9}(index);clear EO_file

%Compute an average altitude above MSL for each VL file
EO.strips(:,1) = unique(EO.vl); %Number of Trimble strips (or individual vl files) 
for mm = 1:length(EO.strips)
    ind_vl = find(EO.strips(mm)==EO.vl);
    EO.alt_avg(ind_vl,1) = mean(EO.alt(ind_vl))-msl;
end

%Compute average cross track resolution for a given camera for each strip/vl file and then compute the nearest 5cm export size for each strip/vl file
if strcmp(cam,'jai')
    dw_x = 18.13; fl = 14; pix_x = 3296; fov_x = 2*atan(dw_x/(2*fl))*180/pi; swath_x = 2*tan(fov_x/2*pi/180)*EO.alt_avg; EO.res_x(:,1) = swath_x/pix_x;
elseif strcmp(cam,'flare')
    dw_x = 22.53; fl = 14; pix_x = 4096; fov_x = 2*atan(dw_x/(2*fl))*180/pi; swath_x = 2*tan(fov_x/2*pi/180)*EO.alt_avg; EO.res_x(:,1) = swath_x/pix_x;
elseif strcmp(cam,'sc6000')
    dw_x = 16; fl = 13; pix_x = 640; fov_x = 2*atan(dw_x/(2*fl))*180/pi; swath_x = 2*tan(fov_x/2*pi/180)*EO.alt_avg; EO.res_x(:,1) = swath_x/pix_x;
elseif strcmp(cam,'sc6700')
    dw_x = 9.6; fl = 13; pix_x = 640; fov_x = 2*atan(dw_x/(2*fl))*180/pi; swath_x = 2*tan(fov_x/2*pi/180)*EO.alt_avg; EO.res_x(:,1) = swath_x/pix_x;
elseif strcmp(cam,'d810')
    dw_x = 35.9; fl = 28; pix_x = 7360; fov_x = 2*atan(dw_x/(2*fl))*180/pi; swath_x = 2*tan(fov_x/2*pi/180)*EO.alt_avg; EO.res_x(:,1) = swath_x/pix_x;
elseif strcmp(cam,'doppvis') %d850 camera on its side so y values become x values!
    dw_x = 24; fl = 20; pix_x = 5504; fov_x = 2*atan(dw_x/(2*fl))*180/pi; swath_x = 2*tan(fov_x/2*pi/180)*EO.alt_avg; EO.res_x(:,1) = swath_x/pix_x;
end
EO.res_cm = ceil(EO.res_x * 100 / 5) * 5+10;
% EO.res_cm(1:end)
EO.res_u = unique(EO.res_cm); clear swath_x

%Write the project summary into a formatted text file
fid = fopen([dirout_prj_xml proj_str(1:end-1) '_Summary.txt'],'w');
fprintf(fid,'-----------------------------------------------------------\n');
fprintf(fid,'| VL File |  UTC Start Time  | Frames | Avg Alt | Exp Res |\n');
fprintf(fid,'-----------------------------------------------------------\n');

%Write vl info in the summary file
for jj=1:length(EO.strips)
    ind = find(EO.strips(jj) == EO.vl);
    time = datenum(proj_str(1:8),'yyyymmdd')+rem(EO.time(ind(1)),24*60*60)/(24*60*60)-18/(24*60*60);
    vl_start = datestr(time,'yyyy/mm/dd HH:MM');
    fprintf(fid,'| %7d | %s | %6d | %5.0fft | %5dcm |\n',EO.strips(jj), vl_start, length(ind), round(EO.alt_avg((ind(1)))*.0328084)*100, EO.res_cm(ind(1)));
end

%Write total frames and total frames per resolution in the summary file
fprintf(fid,'-----------------------------------------------------------\n\n');
fprintf(fid,'Total Flight Images: %d\n',length(EO.frame));
for jj=1:length(EO.res_u)
    ind = find(EO.res_u(jj) == EO.res_cm);
    fprintf(fid,'Total Flight %dcm Images: %d\n',EO.res_u(jj),length(ind));
end
fclose(fid);

%% 2nd NODE: Applications Master PRJ File Generation

%Store the template PRJ file in memory and find all the lines that need to be edited for automation
fid = fopen([dirin_prj 'Template\' temp_prj.name],'r');
ind1 = 0; ind2 = 0; ind3 = 0;

while ~feof(fid)
    ind1 = ind1 + 1; 
    template_prj{ind1} = fgetl(fid);    
    %indices for modifying/adding new text
    if contains(template_prj{ind1},'$PROJECT_NAME') || contains(template_prj{ind1},'$STARTING_DATE') || contains(template_prj{ind1},'$LAST_CHANGE') || contains(template_prj{ind1},'$STRIPS') || contains(template_prj{ind1},'$END_COORDINATE_SYSTEM') || strcmp(template_prj{ind1},'$NAVIGATION')
        ind2 = ind2 + 1; 
        index_ad(ind2) = ind1;
    end    
    %indices for copying information from the template file
    if contains(template_prj{ind1},'$ID') || contains(template_prj{ind1},'$FOCAL_LENGTH') || contains(template_prj{ind1},'$INS_BORESIGHT_ALIGNMENT')
        ind3 = ind3 + 1; 
        index_cp(ind3) = ind1;
    end
end
fclose(fid);

%Retrieve data used for some of the subsequent lines
col_ind = strfind(template_prj{index_cp(1)},':'); %Find colon before camera ID
crs  = template_prj{index_cp(1)}(col_ind+2:end); %Retrieve camera ID name
col_ind = strfind(template_prj{index_cp(2)},':'); %Find colon before camera ID
cam_id = template_prj{index_cp(2)}(col_ind+2:end); %Retrieve camera ID name
col_ind = strfind(template_prj{index_cp(4)},':'); %Find colon before camera focal length
cam_fl = template_prj{index_cp(4)}(col_ind+2:end); %Retrieve camera focal length
col_ind = strfind(template_prj{index_cp(5)},':'); %Find colon before boresight angles
bore_ang = split(template_prj{index_cp(5)}(col_ind+2:end)); %Retrieve camera boresight angles
col_ind = strfind(template_prj{index_cp(6)},':'); %Find colon before navigation ID
nav_id = template_prj{index_cp(6)}(col_ind+2:end); %Retrieve navigation ID

timestamp1 = datestr(now,'ddd mmm dd HH:MM:SS yyyy'); %Current timestamp to add to file
timestamp2 = datestr(now,'HH:MM:SS dd/mm/yyyy'); %Current timestamp to add to file (slightly different format)

%Replace or add lines to PRJ file as needed
template_prj{index_ad(1)} = [template_prj{index_ad(1)} proj_str(1:end-1)]; %$PROJECT_NAME
template_prj{index_ad(2)} = [template_prj{index_ad(2)} timestamp1]; %$STARTING_DATE
template_prj{index_ad(3)} = [template_prj{index_ad(3)} timestamp1]; %$LAST_CHANGE

%Boresight Exterior Orientation Matrix
bo_matrix = [cos(str2num(bore_ang{end-1})*pi/180)*cos(str2num(bore_ang{end})*pi/180),-cos(str2num(bore_ang{end-1})*pi/180)*sin(str2num(bore_ang{end})*pi/180),sin(str2num(bore_ang{end-1})*pi/180);...
    cos(str2num(bore_ang{end-2})*pi/180)*sin(str2num(bore_ang{end})*pi/180)+sin(str2num(bore_ang{end-2})*pi/180)*sin(str2num(bore_ang{end-1})*pi/180)*cos(str2num(bore_ang{end})*pi/180),cos(str2num(bore_ang{end-2})*pi/180)*cos(str2num(bore_ang{end})*pi/180)-sin(str2num(bore_ang{end-2})*pi/180)*sin(str2num(bore_ang{end-1})*pi/180)*sin(str2num(bore_ang{end})*pi/180),-sin(str2num(bore_ang{end-2})*pi/180)*cos(str2num(bore_ang{end-1})*pi/180);...
    sin(str2num(bore_ang{end-2})*pi/180)*sin(str2num(bore_ang{end})*pi/180)-sin(str2num(bore_ang{end-1})*pi/180)*cos(str2num(bore_ang{end-2})*pi/180)*cos(str2num(bore_ang{end})*pi/180),sin(str2num(bore_ang{end-2})*pi/180)*cos(str2num(bore_ang{end})*pi/180)+cos(str2num(bore_ang{end-2})*pi/180)*sin(str2num(bore_ang{end-1})*pi/180)*sin(str2num(bore_ang{end})*pi/180),cos(str2num(bore_ang{end-2})*pi/180)*cos(str2num(bore_ang{end-1})*pi/180)];

%Add $STRIPS PRJ line information. Suppressed for now, may need later.
% for jj=1:length(EO.strips)   
%     index = find(EO.vl == EO.strips(jj));
%     str1 = sprintf('    %-3d   ElementPhoto   0.00 { ',EO.strips(jj)); 
%     str2 = '';
%     for kk=1:length(index)
%         str2 = [str2 sprintf('%03d%04d ',EO.vl(index(kk)),EO.frame(index(kk)))];
%     end    
%     str3 = '}';
%     strip{jj,:} = [str1 str2 str3];  
% end


%Add $PHOTO PRJ line information
for mm=1:length(EO.frame)    
    %Enter basic information
    photo{mm,1,:} = '$PHOTO'; 
    photo{mm,2,:} = sprintf('  $PHOTO_NUM : %03d%04d',EO.vl(mm),EO.frame(mm));
    photo{mm,3,:} = sprintf('  $PHOTO_FILE : %s',[dirin_img filelist(mm).name]);
    photo{mm,4,:} = sprintf('  $CAMERA_ID : %s',cam_id);
    photo{mm,5,:} = sprintf('  $TERRAIN_HEIGHT : %f',msl);
    photo{mm,6,:} = '  $IO_STAT : manual';
    photo{mm,7,:} = '  $ORI_STAT : Initial';
    photo{mm,8,:} = sprintf('  $EXT_ORI : %s',timestamp2);
    photo{mm,9,:} = sprintf('      %10.5f   %13.5f   %13.5f   %13.5f',str2num(cam_fl),EO.east(mm),EO.north(mm),EO.alt(mm));    
    %Calculate exterior orientation from boresight angles and EO file rounded values and write to cell array
    eo_matrix = [cos(round(EO.phi(mm),5)*pi/180)*cos(round(EO.kappa(mm),5)*pi/180),-cos(round(EO.phi(mm),5)*pi/180)*sin(round(EO.kappa(mm),5)*pi/180),sin(round(EO.phi(mm),5)*pi/180);...
        cos(round(EO.omega(mm),5)*pi/180)*sin(round(EO.kappa(mm),5)*pi/180)+sin(round(EO.omega(mm),5)*pi/180)*sin(round(EO.phi(mm),5)*pi/180)*cos(round(EO.kappa(mm),5)*pi/180),cos(round(EO.omega(mm),5)*pi/180)*cos(round(EO.kappa(mm),5)*pi/180)-sin(round(EO.omega(mm),5)*pi/180)*sin(round(EO.phi(mm),5)*pi/180)*sin(round(EO.kappa(mm),5)*pi/180),-sin(round(EO.omega(mm),5)*pi/180)*cos(round(EO.phi(mm),5)*pi/180);...
        sin(round(EO.omega(mm),5)*pi/180)*sin(round(EO.kappa(mm),5)*pi/180)-sin(round(EO.phi(mm),5)*pi/180)*cos(round(EO.omega(mm),5)*pi/180)*cos(round(EO.kappa(mm),5)*pi/180),sin(round(EO.omega(mm),5)*pi/180)*cos(round(EO.kappa(mm),5)*pi/180)+cos(round(EO.omega(mm),5)*pi/180)*sin(round(EO.phi(mm),5)*pi/180)*sin(round(EO.kappa(mm),5)*pi/180),cos(round(EO.omega(mm),5)*pi/180)*cos(round(EO.phi(mm),5)*pi/180)];   
    eobo_matrix = (eo_matrix*bo_matrix)';    
    photo{mm,10,:} = sprintf('      %15.12f      %15.12f      %15.12f',eobo_matrix(1,1),eobo_matrix(1,2),eobo_matrix(1,3));
    photo{mm,11,:} = sprintf('      %15.12f      %15.12f      %15.12f',eobo_matrix(2,1),eobo_matrix(2,2),eobo_matrix(2,3));
    photo{mm,12,:} = sprintf('      %15.12f      %15.12f      %15.12f',eobo_matrix(3,1),eobo_matrix(3,2),eobo_matrix(3,3));
    photo{mm,13,:} = '$END';  
end  

%Add $NAVIGATION PRJ line information
for qq=1:length(EO.frame)    
    nav{qq,1,:} = sprintf('  $PHOTO_NUM : %03d%04d',EO.vl(qq),EO.frame(qq));
    nav{qq,2,:} = sprintf('    $NAVIGATION_PAR : %s',nav_id);
    nav{qq,3,:} = sprintf('    $GPS :     %11.3f     %11.3f     %10.3f',EO.east(qq),EO.north(qq),EO.alt(qq));
    nav{qq,4,:} = sprintf('    $INS :     %10.5f     %10.5f     %10.5f',EO.omega(qq),EO.phi(qq),EO.kappa(qq));
    nav{qq,5,:} = sprintf('    $COORDINATE_SYSTEM_ID : %s',crs);      
end 

%Write new file using the template with added strip, photo, and nav data for each output resolution
for yy=1:length(EO.res_u)    
    ind_res = find(EO.res_cm == EO.res_u(yy)); %find indices of photos at prescribed resolution
    fid = fopen([dirout_prj_xml proj_str num2str(EO.res_u(yy)) 'cm.prj'],'w');   
    for aa=1:index_ad(4)
        fprintf(fid,'%s\n',template_prj{aa});
    end
    %Suppress strips for now, may use in the future
    % for aa=1:size(strip,1)
    %     fprintf(fid,'%s\n',strip{aa});
    % end
    for aa=index_ad(4)+1:index_ad(5)
        fprintf(fid,'%s\n',template_prj{aa});
    end    
    for tt=1:length(ind_res)
        for aa=1:size(photo,2)
            fprintf(fid,'%s\n',photo{ind_res(tt),aa,:});
        end
    end
    for aa=index_ad(5)+1:index_ad(6)
        fprintf(fid,'%s\n',template_prj{aa});
    end
    for tt=1:length(ind_res)
        for aa=1:size(nav,2)
            fprintf(fid,'%s\n',nav{ind_res(tt),aa,:});
        end
    end
    for aa=index_ad(6)+1:length(template_prj)
        fprintf(fid,'%s\n',template_prj{aa});
    end
    fclose(fid);
end

%% 3rd NODE: Orthomaster XML File Generation

%Store the template XML file in memory and find all the lines that need to be edited for automation
fid = fopen([dirin_prj 'Template\' temp_xml.name],'r');
ind1 = 0; ind2 = 0;
while ~feof(fid)
    ind1 = ind1 + 1; 
    template_xml{ind1} = fgetl(fid);    
    %indices for modifying/adding new text
    if contains(template_xml{ind1},'dehPlane0') || contains(template_xml{ind1},'resolution') || contains(template_xml{ind1},'outdirPath') || contains(template_xml{ind1},'nameMask')
        ind2 = ind2 + 1; 
        index_mod(ind2) = ind1;
    end    
end
fclose(fid);

%Retrieve data used for some of the subsequent lines
dem_ind = strfind(template_xml{index_mod(1)},'"'); %Find first quote to insert dem value
template_xml{index_mod(1)} = sprintf('%s%.2f%s',template_xml{index_mod(1)}(1:dem_ind),msl,template_xml{index_mod(1)}(dem_ind+1:end));
dir_ind = strfind(template_xml{index_mod(3)},';'); %Find first semicolon to insert outpu directory
template_xml{index_mod(3)} = sprintf('%s%s%s',template_xml{index_mod(3)}(1:dir_ind),dirout_ortho(1:end-1),template_xml{index_mod(3)}(dir_ind+1:end));

%write xml file for each unique 5cm-rounded resolution output
for yy=1:length(EO.res_u)
    res_ind = strfind(template_xml{index_mod(2)},'"'); %Find first quote to insert resolution value
    mask_ind = strfind(template_xml{index_mod(4)},'"'); %Find first quote to insert output name mask
    template_xml{index_mod(2)} = sprintf('%s%.2f%s',template_xml{index_mod(2)}(1:res_ind(1)),EO.res_u(yy)/100,template_xml{index_mod(2)}(res_ind(2):end));
    template_xml{index_mod(4)} = sprintf('%s%s%dcm_&lt;PHOTO>%s',template_xml{index_mod(4)}(1:mask_ind(1)),proj_str,EO.res_u(yy),template_xml{index_mod(4)}(mask_ind(2):end));   
    fid = fopen([dirout_prj_xml proj_str num2str(EO.res_u(yy)) 'cm.xml'],'w');
    for aa=1:length(template_xml)
        fprintf(fid,'%s\n',template_xml{aa});
    end
    fclose(fid);
end

%% 4th NODE: Batch File Generation

%Write batch file to run all of the PRJ/XML files for an individual flight (if multiple resolutions exist). Right click and run as Administrator
%Change to the directory of the Orthomaster module and run for each set of PRJ/XML resolution files
trimble_dir = 'C:\Program Files\Trimble\Trimble Photogrammetry 11.0\bin'; 
prj_files = dir([dirout_prj_xml '*cm.prj']);
fid = fopen([dirout_prj_xml proj_str(1:end-1) '.bat'],'w');
fprintf(fid,'cd /d %s\n',trimble_dir);

%Write the total number of images and the total number for each resolution to echo during georectification
fprintf(fid,'echo %d Total Images; ',length(EO.frame));
for kk=1:length(EO.res_u)
    ind = find(EO.res_u(kk) == EO.res_cm);
    fprintf(fid,'%d %dcm Images; ',length(ind),EO.res_u(kk));
end
fprintf(fid,'\n');

%Write the Trimble batch commands for each resolution
for kk=1:length(prj_files)
    fprintf(fid,'orthomaster.exe -batch -prj %s%s\n',dirout_prj_xml,prj_files(kk).name);
end

fprintf(fid,'break> %s%s\n',dirV,'blank.txt');
fclose(fid);



