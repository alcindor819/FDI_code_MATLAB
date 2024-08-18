function user_call(data_name,distribution,P,FDR,res_path,X,Y,res_name,chr_n,numWorkers)
% data_name = 'D:/wyzwork/Basic_info/IH02_data/';
% distribution = 'Beta';
% P = 0.01;
% FDR = 0.05;
% res_path = 'D:/wyzwork/精简版代码/FDI/Dispersed_region_identify/call/';
% res_name = 'IH02_test';
% X = 0.5;
% Y = 10;
% chr_n = 1:22;
% numWorkers = 2;%核数取2，约1.5小时运行完
parpool('local', numWorkers);
parfor i = chr_n
    Dispersed_region_identify_chr(data_name,distribution,P,FDR,res_path,X,Y,res_name,i);
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