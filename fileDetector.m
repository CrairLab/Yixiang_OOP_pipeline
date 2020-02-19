function fileDetector()
%{
    Generate txt files based on file detected. Generate files.txt as a list
    of names of .tif movies. Generate Spike2files.txt as a list of names of 
    .mat spike2 data. Generat baphyfiles.txt as a list of names of baphy
    files. 
%}
   
    tmp_list = dir;
    screen_types = {'.mat','.m','.tif'};
    disp(['Working Folder: ' pwd])
    tifID = fopen('files.txt','w+');
    spike2ID = fopen('Spike2files.txt','w+');
    baphyID = fopen('baphyfiles.txt','w+');
    
    for i = 1 : length(tmp_list)
        
        curName = tmp_list(i).name;
        
        if length(curName) < 2
            continue
        end
        
        last2chars_name = curName(end-1:end);
        count = 0;
        flag = false;
        
        while (count < length(screen_types))&&(~flag)
            count = count + 1;
            curr_type = screen_types{count};
            last2chars_type = curr_type(end-1:end);
            
            if last2chars_name == last2chars_type
                flag = true;
                switch curr_type
                    case '.mat'
                        fprintf(spike2ID,'%s\r\n',curName);
                    case '.m'
                        fprintf(baphyID,'%s\r\n',curName);
                    case '.tif'
                        if isempty(strfind(curName,'@00'))
                            fprintf(tifID,'%s\r\n',curName);
                        end
                end
            end
                    
        end
    end
    
    fclose('all');
end