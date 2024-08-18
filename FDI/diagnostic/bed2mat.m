function bed2mat(data_name,f1,f3,d3)%将bed文件转换为22个chr.mat

%举例：bed222( 'HL1.Sample1.arrayexpress.bed','/home/public_data/cfDNA/chehuiwen/unzip_pancancer/','/home/public_data/cfDNA/chehuiwen/mat_pancancer/','Sample 1');
%data_name是bed文件的名字如上第一个，f1是bedpath，bed文件存放的文件夹，f3是输出的文件夹名字
%该函数会在f3文件夹里生成每个样本的文件夹，d3是生成样本文件夹的名字（自定义）
file_name1 = f1;%
file_name2 = strcat(file_name1,data_name);
file_name3 = f3;%outpath
fileID = fopen(file_name2);
C = textscan(fileID,'%s %d %d %s %d %c', 'commentStyle', '#');  %%import the bed file of the reads
fclose(fileID);
%这里注意，读不同的Bed，每一列参数可能不一样，上面textscan可能需要改参数，可以自己查对应参数，第五列是匹配质量，一般保留30以上的，没有就不管
c11 = C{1,1};


file_path3 = strcat(file_name3,d3);
mkdir (file_path3);
file_path4 = strcat(file_path3,'/');
f5 = strcat(file_path4,'data_n/');
mkdir(f5);


ic = 1;
read1 = zeros(length(c11),2);
i1 = 1;

dqpath = pwd;
positions = strfind(dqpath, '/');
lastPosition = positions(end);
lastpath = dqpath(1:lastPosition);
bipath = strcat(lastpath,'Basic_info/');
info_path = strcat(bipath,'chr_info.mat');
Sinfo = load(info_path);
chrinfo = Sinfo.chrinfo;
dfpath = strcat(bipath,'dark_flag/');

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
        Sdl = load(read_dfath);
        dark_flag = Sdl.darkflag;
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


