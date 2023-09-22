function cleanEOfile(EOpath,EOdir,tracks,D)
for i=1:length(tracks)
    if (~isempty(tracks(i).Indices))
    ind = strfind(D(1).name,'_');
    tempStr=[EOdir D(1).name(1:ind(end)) 'track_' num2str(i) '_' D(1).name(1+ind(end):end)];
    
    fidN = fopen(tempStr,'w');
    fid=fopen(EOpath);
    tline = fgetl(fid);
    fprintf(fidN,'%s\n',tline);
    Counter=1;
    while ~contains(tline,'(sec)')
        tline = fgetl(fid);
        Counter=Counter+1;
        fprintf(fidN,'%s\n',tline);
    end
    
    
    
    file=fileread(EOpath);
    lines = regexp(file, '\n', 'split');
    
    for j=tracks(i).Indices(1):tracks(i).Indices(2)
        fprintf(fidN,'%s\n',lines{1,j+Counter});
    end
    
    fclose(fid);
    fclose(fidN);
    end
end