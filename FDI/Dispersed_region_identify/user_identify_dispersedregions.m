% bed_folder = 'D:\wyzwork\精简版代码\FDI\Dispersed_region_identify\bed_folder\';
% bed_name = 'SRR16574631';
% mat_path = 'D:\wyzwork\精简版代码\FDI\Dispersed_region_identify\mat\';
% distribution = 'Beta';
% P = 0.01;
% FDR = 0.05;
% res_path = 'D:\wyzwork\精简版代码\FDI\Dispersed_region_identify\call_new\';
% X = 0.5;
% Y = 10;
% chr_n = 1:22;
% numWorkers = 2;

function user_identify_dispersedregions(bed_folder,mat_path,distribution,P,FDR,res_path,X,Y,chr_n,numWorkers)
bed_name_bed = strcat(bed_name,'.bed');
bed2mat(bed_name_bed,bed_folder,mat_path,bed_name);
res_name =bed_name;
data_name = strcat(mat_path,bed_name);
data_name = strcat(data_name,'\');
user_call(data_name,distribution,P,FDR,res_path,X,Y,res_name,chr_n,numWorkers);
res_path2 = strcat(res_path,res_name);
res_path2 = strcat(res_path2,'\');
Dispersed_region_analysis(res_path2);
end