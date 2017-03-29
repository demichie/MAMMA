function read_bak(pathname,filename)

fid = fopen(strcat(pathname,filename));
fid2= fopen('tmp_file','w');

for i=1:200,
    data_line = fgetl(fid);
    if (data_line==-1) 
        break
    else
        data_line = strrep(data_line,'=',' = ');
        fprintf(fid2,' %s \n', data_line);
    end
end

fclose(fid2);
fclose(fid);

fid = fopen('tmp_file','r');

tline = textscan(fid, '%s');

len = length(tline{1}); 

for i = 1:len,
   
    b = strrep(tline{1}(i),',','');
    b = strrep(b,'"','');
    b = strrep(b,'''','');
    

    if ( isempty(b{1}) )
        
        len = len - 1;
        
        tline{1}(i:len) = tline{1}(i+1:len+1);
        
    else
        
        tline{1}(i) = b;
        
    end
            
end


strings = tline{1}(1:len);

set_variables(strings);

fclose(fid);

function set_variables(strings)

len = length(strings);

for i= 1:len,
    
    testo = strings{i};
    C = strcmp(testo,'=');
    
    if C == 1

        strings{i-1} = strrep(strings{i-1},'%','_');
        
        if ( isnan( str2double(strings{i+1}) ) )
            
            assignin('base',strings{i-1},strings{i+1});

        else

            assignin('base',strings{i-1},str2double(strings{i+1}));

        end

    end
end



