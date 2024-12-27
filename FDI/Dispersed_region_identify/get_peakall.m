function get_peakall(data_name,distribution,P,FDR,res_path,X,Y,res_name,chr_n,numWorkers,rds,csz,csy)

parpool('local', numWorkers);
parfor i = chr_n
    Dispersed_region_identify_chr(data_name,distribution,P,FDR,res_path,X,Y,res_name,i,rds,csz,csy)
end
peak_all = zeros(1,8);
for i = chr_n
    file_name  = strcat(res_path,res_name);
    file_name = strcat(file_name,'/result_n/peak_');
    file_name = strcat(file_name,num2str(i));
    file_name = strcat(file_name,'.mat');
    load(file_name);
    peak_all = [peak_all;peak_a];
end
peak_all(1,:) = [];
save_path = strcat(res_path,res_name);
save_path = strcat(save_path,'/result_n/');
peak_a = peak_all;
save([save_path,'peak_all'],'peak_a');
end