function bed2mat(data_name,f1,f3,d3)%Converting bed files to 22 chr.mat


file_name1 = f1;%bed folder path
file_name2 = strcat(file_name1,data_name);
file_name3 = f3;%outpath
fileID = fopen(file_name2);
C = textscan(fileID,'%s %d %d %s %d %c', 'commentStyle', '#');  %%import the bed file of the reads
fclose(fileID);
%Note here, read different bed, each column parameters may not be the same, the above textscan may need to change the parameters,
%you can check the corresponding parameters, the fifth column is the quality of the match, generally retained more than 30, no matter if there is no
c11 = C{1,1};
chrinfo = cell(22,2);
file_path3 = strcat(file_name3,d3);
mkdir (file_path3);
file_path4 = strcat(file_path3,'/');
f5 = strcat(file_path4,'data_n/');
mkdir(f5);
ic = 1;
read1 = zeros(length(c11),2);
i1 = 1;
dqpath = pwd;
positions = strfind(dqpath, '\');
lastPosition = positions(end);
lastpath = dqpath(1:lastPosition);
bipath = strcat(lastpath,'Basic_info\');
info_path = strcat(bipath,'chr_info.mat');
load(info_path);
dfpath = strcat(bipath,'dark_flag\');

for izong = 1:length(c11)
    name = cell2mat(C{1,1}(izong));
    name = strrep(name,'chr','');
    if strcmp(name,chrinfo{ic,1})
        if C{1,5}(izong)>29
            read1(i1,1) = C{1,2}(izong)+1;
            read1(i1,2) = C{1,3}(izong)-C{1,2}(izong)+1;
            i1 = i1+1;
        end
    else
        region_len = cell2mat(chrinfo(ic,2));
        i1 = 1;
        ic1 = num2str(ic);
        read_dfath = strcat(dfpath,'chr');
        read_dfath = strcat(read_dfath,ic1);
        read_dfath = strcat(read_dfath,'.mat');
        load(read_dfath);
        clear read;
        read = read1;
        read(~any(read,2),:)=[];
        chrn1 = cell2mat(chrinfo(ic,1));
        chrn2 = strcat('chr',chrn1);
        save([f5,chrn2],'read','region_len','dark_flag');
        clear read1;
        
        
        namee = str2double(name);
        na = isnan(namee);
        if na==0
            ic = namee;
        else
            break;
        end
        if C{1,5}(izong)>29
            read1(i1,1) = C{1,2}(izong)+1;
            read1(i1,2) = C{1,3}(izong)-C{1,2}(izong)+1;
            i1 = i1+1;
        end
    end
end
end


