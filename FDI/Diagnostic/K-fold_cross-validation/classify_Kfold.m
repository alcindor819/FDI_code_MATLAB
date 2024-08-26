

function classify_Kfold(sample_info_path,name_col,pos_neg_col,fold,mat_folder,X,Y,readlen,buchang,gp,fdr)

[~,~,info] = xlsread(sample_info_path);
sample_name = info(:,name_col);
sample_category = cell2mat(info(:,pos_neg_col));

c1 = sample_name((sample_category==0),:);
h1 = sample_name((sample_category==1),:);
dataname = [c1;h1];
dqpath=pwd;%当前文件夹路径
dqpath = strcat(dqpath,'/');
rng('default');
N = length(dataname);%RNG DEFULT
n1 = zeros(N,1);
for ii = 1:size(dataname,1)
    if ismember(dataname(ii),h1)
        n1(ii) = 1;
    else
        n1(ii) = 0;
    end
end
Indices = crossvalind('Kfold', n1, fold);   %%3交叉验证，将样本分成3份
cfdnapath = mat_folder;

positions = strfind(dqpath, '/');
lastPosition = positions(end-2);
lastpath = dqpath(1:lastPosition);
bipath = strcat(lastpath,'Basic_info/');
map_path = strcat(bipath,'mappability/');

parfor j =1:fold
    Kfold_FDI(j,c1,h1,Indices,dqpath,cfdnapath,X,Y,readlen,buchang,map_path,gp,fdr);
end
for i =1:fold
    i1 = num2str(i);
    ddpath = strcat(dqpath,i1);
    ddpath = strcat(ddpath,'/');
    label  = ismember(dataname,h1);
    
    for i9 =1:length(dataname)
        ln = strcat(ddpath,cell2mat(dataname(i9)));
        ln1 = strcat(ln,'/feature.mat');
        
        load(ln1);
        feature_all(:,i9) = feature3;
        ln2 = strcat(ln,'/feature2000.mat');
        
        load(ln2);
        feature_all2000(:,i9) = feature3;
    end
    data = feature_all;
    clear feature_all;
    data(isnan(data)) = 0;
    data(:,all(data==0,1))=[];
    
    data2000 = feature_all2000;
    clear feature_all2000;
    data2000(isnan(data2000)) = 0;
    data2000(:,all(data2000==0,1))=[];
    
    
    
    for i33 = 1:size(data,2)
        data(:,i33) = zscore(data(:,i33));
    end
    
    
    
    
    for i33 = 1:size(data2000,2)
        data2000(:,i33) = zscore(data2000(:,i33));
    end
    
    j = i;
    train_data=data(:,Indices~=j)';     %%训练集数据
    train_data(isnan(train_data)) = 0;
    
    train_label=double(label(Indices~=j,1));  %%训练集标签
    
    test_data=data(:,Indices==j)';     %%测试集数据
    test_data(isnan(test_data)) = 0;
    
    test_label=double(label(Indices==j,1));  %%测试集标签
    
    %
    train_data2000=data2000(:,Indices~=j)';     %%训练集数据
    train_data2000(isnan(train_data2000)) = 0;
    
    train_label=double(label(Indices~=j,1));  %%训练集标签
    
    test_data2000=data2000(:,Indices==j)';     %%测试集数据
    test_data2000(isnan(test_data2000)) = 0;
    
    test_label=double(label(Indices==j,1));  %%测试集标签
    
    
    
    
    nb_svm=fitcsvm(train_data,train_label);    %%SVM的训练模型
    [res_svm, score] = predict(nb_svm,test_data);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分
    
    
    [X,Y,~,auc]=perfcurve(test_label,score(:,2),'1'); %计算分类器的AUC
    intervals= linspace(0, 1, 101);
    val=data_point_estimate(X,Y,intervals);
    if j==1
        mean_curve= val/fold_number;
    else
        mean_curve= mean_curve+ val/fold_number;
    end
    
    
    
    %
    nb_svm2=fitcsvm(train_data2000,train_label);    %%SVM的训练模型
    [res_svm, score] = predict(nb_svm2,test_data2000);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分
    
    
    [X2,Y2,~,auc]=perfcurve(test_label,score(:,2),'1'); %计算分类器的AUC
    intervals= linspace(0, 1, 101);
    val2=data_point_estimate(X2,Y2,intervals);
    if j==1
        mean_curve2= val2/fold_number;
    else
        mean_curve2= mean_curve2+ val2/fold_number;
    end
    
    
    
end
X=[0 intervals]';
Y=[0;mean_curve];

X2=[0 intervals]';
Y2=[0;mean_curve2];
d2 = dqpath;
save([d2,'roc'] ,'X','Y');
save([d2,'roc2'] ,'X2','Y2');

fold = 10;
sum = 0;
sum2000 = 0;
dqpath = pwd;
dqpath = strcat(dqpath,'/');
for i =1:fold
    
    i1 = num2str(i);
    dqpath = strcat(dqpath,i1);
    dqpath = strcat(dqpath,'/');
    i11 = strcat(i1,'_2000.mat');
    i2 = strcat(i1,'.mat');
    i111 = strcat(dqpath,i11);
    i22 = strcat(dqpath,i2);
    load(i111);
    sum2000 = sum2000+cell2mat(res(2,1));
    load(i22);
    sum = sum+cell2mat(res(2,1));
    dqpath = pwd;
    dqpath = strcat(dqpath,'/');
end
auc2000 = sum2000/fold;
auc = sum/fold;
auc = [auc auc2000];
save([dqpath,'auc'] ,'auc');
end