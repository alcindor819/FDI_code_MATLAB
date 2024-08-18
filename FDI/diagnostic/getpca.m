function [auc100,auc100_2000] = getpca(j,c1,h1,Indices,dqpath)%这是交叉验证的


fname = num2str(j);%第几折，整个文件夹，dqpath代表这折的文件夹


dqpath = strcat(dqpath,fname);
dqpath = strcat(dqpath,'/');

dataname = [c1;h1];%总数据样本名字
train_name = dataname(Indices~=j);%综合训练集名字
%test_name = dataname(Indices==j);%综合测试集名字


label  = ismember(dataname,h1);%癌症即阳性为1
for i9 =1:length(dataname)
    ln = strcat(dqpath,cell2mat(dataname(i9)));
    ln = strcat(ln,'/feature.mat');
    load(ln);
    feature_all(:,i9) = feature3;
end
for i9 =1:length(dataname)
    ln = strcat(dqpath,cell2mat(dataname(i9)));
    ln = strcat(ln,'/feature2000.mat');
    load(ln);
    feature_all2000(:,i9) = feature3;
end
data = feature_all;
data2000 = feature_all2000;
data(isnan(data)) = 0;
data(:,all(data==0,1))=[];
data2000(isnan(data2000)) = 0;
data2000(:,all(data2000==0,1))=[];
for i33 = 1:size(data,2)
    data(:,i33) = zscore(data(:,i33));
end
for i33 = 1:size(data2000,2)
    data2000(:,i33) = zscore(data2000(:,i33));
end
%开始照搬分类器原始peak
train_data=data(:,Indices~=j)';     %%训练集数据，一行一样本，一列一区间
train_data2000=data2000(:,Indices~=j)';
train_data(isnan(train_data)) = 0;
train_data2000(isnan(train_data2000)) = 0;
train_label=double(label(Indices~=j,1));  %%训练集标签，一行一样本

test_data=data(:,Indices==j)';     %%测试集数据
test_data(isnan(test_data)) = 0;
test_data2000=data2000(:,Indices==j)';     %%测试集数据
test_data2000(isnan(test_data2000)) = 0;
test_label=double(label(Indices==j,1));  %%测试集标签

x1 = train_data;
x2 = train_data2000;
y1 = test_data;
y2 = test_data2000;

auc100 = zeros(100,2);
auc100_2000 = zeros(100,2);
% 对特征矩阵进行PCA降维
[coeff, ~, latent] = pca(x1);
[coeff2, ~, latent2] = pca(x2);
% 计算每个主成分的方差贡献率
variances = cumsum(latent) / sum(latent);
variances2 = cumsum(latent2) / sum(latent2);
% 选择方差贡献率较高的主成分
for i = 1:100
    variance_threshold = i*0.01;
    selected_components = find(variances >= variance_threshold, 1);
    
    % 将特征矩阵投影到选定的主成分上
    x22 = x1 * coeff(:, 1:selected_components);
    y22 = y1 * coeff(:, 1:selected_components);
    
    
    
    
    
    nb_svm2=fitcsvm(x22,train_label);
    [res_svm2, score2] = predict(nb_svm2,y22);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分
    
    [X,Y,~,auc]=perfcurve(test_label,score2(:,2),'1'); %计算分类器的AUC
    auc100(i,1) = auc;
    auc100(i,2) = size(x22,2);
end
for i = 1:100
    variance_threshold = i*0.01;
    selected_components = find(variances2 >= variance_threshold, 1);
    
    % 将特征矩阵投影到选定的主成分上
    x22 = x2 * coeff2(:, 1:selected_components);
    y22 = y2 * coeff2(:, 1:selected_components);
    
    
    
    
    
    nb_svm2=fitcsvm(x22,train_label);
    [res_svm2, score2] = predict(nb_svm2,y22);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分
    
    [X,Y,~,auc]=perfcurve(test_label,score2(:,2),'1'); %计算分类器的AUC
    auc100_2000(i,1) = auc;
    auc100_2000(i,2) = size(x22,2);
end
save([dqpath,'auc100'],'auc100','auc100_2000');

end