function user_processing(bed_info_path,bed_folder,numWorkers,res_folder,name_col)
[~,~,r] = xlsread(bed_info_path);
parpool('local', numWorkers);
parfor i = 1:size(r,1)
    name = cell2mat(r(i,name_col));
    name_bed = strcat(name,'.bed');
    bed2mat(name_bed,bed_folder,res_folder,name);
end
end
