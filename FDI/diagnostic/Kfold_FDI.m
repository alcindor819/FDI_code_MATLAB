function Kfold_FDI(j,c1,h1,Indices,dqpath,cfdnapath,X,Y,readlen,buchang,map_path,gp,fdr)%优化了合并速度，但是算特征很慢，可以判断peak存在

fname = num2str(j);%第几折，整个文件夹，dqpath代表这折的文件夹
mkdir(fname);
jk_dirname = strcat(fname,'/jk_traindata/');
mkdir(jk_dirname);
hcc_dirname = strcat(fname,'/hcc_traindata/');
mkdir(hcc_dirname);
dqpath = strcat(dqpath,fname);
dqpath = strcat(dqpath,'/');


x = 20; % 假设需要选择 30 行
if x>=floor(0.8*min(size(h1,1),size(c1,1)))
    x=floor(0.8*min(size(h1,1),size(c1,1)));
end
selected_rows = randsample(size(h1, 1), x);
selected_rows2 = randsample(size(c1, 1), x);
% 获取选择的行
selected_h11 = h1(selected_rows, :);
selected_c11 = c1(selected_rows2, :);
h2 = selected_h11;
c2 = selected_c11;



dataname = [c1;h1];%总数据样本名字



    
    for i =1:length(dataname)%生成对应的文件夹，其实就存放以下对应的特征
        lj1 = dqpath;
        lj2 = cell2mat(dataname(i));
        lj1 = strcat(lj1,lj2);
        lj1 = strcat(lj1,'/');
        mkdir(lj1);
    end
    for i =1:22
        ii = num2str(i);
        name = strcat('chr',ii);
        read = zeros(1,2);
        
        save([jk_dirname, name],'read');
    end%生成初始chr
    for i =1:22
        ii = num2str(i);
        name = strcat('chr',ii);
        read = zeros(1,2);
        save([hcc_dirname, name],'read');
    end
    parfor i7 = 1:22%把新的CFDNA加进去
        chrname = strcat('chr',num2str(i7));
        chrname = strcat(chrname,'.mat');
        l3 = strcat(jk_dirname,chrname);
        s = load(l3);
        read_0 = s.read;
        for i1 = 1:x%整合JK训练集的CFDNA
            loadname = cell2mat(c2(i1));
            loadname = strcat(cfdnapath,loadname);
            loadname = strcat(loadname,'/data_n/');
            
            loadname1 = strcat(loadname,chrname );
            
            s1 = load(loadname1);%读取将要合并的数据
            read_0 = [read_0;s1.read];
        end
        read_0(~any(read_0,2),:)=[];
        read_0 = sortrows(read_0);
        read = read_0;
        savename = jk_dirname;
        r_parsave(savename,chrname,read,s1.region_len,s1.dark_flag);
        
    end
    parfor i7 = 1:22%把新的CFDNA加进去
        chrname = strcat('chr',num2str(i7));
        chrname = strcat(chrname,'.mat');
        l3 = strcat(hcc_dirname,chrname);
        s = load(l3);
        read_0 = s.read;
        for i1 = 1:x%整合JK训练集的CFDNA
            loadname = cell2mat(h2(i1));
            loadname = strcat(cfdnapath,loadname);
            loadname = strcat(loadname,'/data_n/');
            
            loadname1 = strcat(loadname,chrname );
            
            s2 = load(loadname1);%读取将要合并的数据
            read_0 = [read_0;s2.read];
        end
        read_0(~any(read_0,2),:)=[];
        read_0 = sortrows(read_0);
        read = read_0;
        savename = hcc_dirname;
        r_parsave(savename,chrname,read,s2.region_len,s2.dark_flag);
        
    end
    clear read;
    clear read0;
    parfor i5 =1:2%并行去call
        if i5==1
            FDI_call(gp,fdr,jk_dirname,jk_dirname,'jkpeak',X,Y,readlen,buchang,map_path);
        else
            FDI_call(gp,fdr,hcc_dirname,hcc_dirname,'hccpeak',X,Y,readlen,buchang,map_path);
        end
    end
    
    
    hebingjkpeak = strcat(jk_dirname,'jkpeak/result_n/peak_all.mat');
    load(hebingjkpeak);
    peak_a1 = peak_a;
    hebinghccpeak = strcat(hcc_dirname,'hccpeak/result_n/peak_all.mat');
    load(hebinghccpeak);
    peak_a = [peak_a1;peak_a];
    
    mm = size(peak_a,1);
    if mm>2
        for i=1:size(peak_a,1)-1
            if peak_a(i+1,2)<peak_a(i,3) && peak_a(i+1,1)==peak_a(i,1)%连接一下互相重叠的
                peak_a(i+1,2)=peak_a(i,2);
                peak_a(i+1,3) = max(peak_a(i+1,3),peak_a(i,3));
                peak_a(i+1,4)=max(peak_a(i+1,4),peak_a(i,4));
                peak_a(i+1,7)=max(peak_a(i+1,7),peak_a(i,7));
                peak_a(i+1,5)=max(peak_a(i+1,5),peak_a(i,5));
                peak_a(i+1,6)=max(peak_a(i+1,6),peak_a(i,6));
                peak_a(i+1,8)=peak_a(i+1,3)-peak_a(i+1,2);
                peak_a(i,:)=0;
            end
        end%连接重叠的
    end
    peak_a(~any(peak_a,2),:)=[];
    peak_a = sortrows(peak_a,1);
    save([dqpath,'peak_all'],'peak_a');
    

parfor ij = 1:length(dataname)
    
    i55 = cell2mat(dataname(ij));%具体人的文件夹名字
    
    feature2000 = getfeature2000(i55,dqpath,cfdnapath,X,Y);%自己写的算特征函数，因为变量在Parfor里不共享且可能出错
    ln2 = strcat(dqpath,i55);
    ln2 = strcat(ln2,'/');
    parsave(feature2000(:,1),ln2,'feature');%自己写的存储特征函数，因为不能用save
    parsave(feature2000(:,2),ln2,'feature2000');
end

%特征都有了，岂不是应该分类了
%先把特征放一起，顺便产生一个标签
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
train_data(isnan(train_data)) = 0;
train_label=double(label(Indices~=j,1));  %%训练集标签，一行一样本

test_data=data(:,Indices==j)';     %%测试集数据
test_data(isnan(test_data)) = 0;
test_label=double(label(Indices==j,1));  %%测试集标签

nb_svm=fitcsvm(train_data,train_label);    %%SVM的训练模型
[res_svm, score] = predict(nb_svm,test_data);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分

[~,~,~,auc]=perfcurve(test_label,score(:,2),'1'); %计算分类器的AUC

auc_svm=auc;
cpre=classperf(test_label,res_svm);  %计算分类器的准确率等指标
Accuracy_svm=cpre.CorrectRate;  %%准确率
sensitivity_svm=cpre.Sensitivity;  %%灵敏度
specificity_svm=cpre.Specificity;   %%特异度
result=[auc_svm Accuracy_svm sensitivity_svm specificity_svm]./1;

title={'AUC','Accuracy','Seneitivity','Specificity'};
res=[title;num2cell(result)];
ressavename = strcat(dqpath,fname);
ressavename = strcat(ressavename,'.mat');
save(ressavename,'res');
%开始照搬分类器   2000peak
train_data2000=data2000(:,Indices~=j)';     %%训练集数据
train_data2000(isnan(train_data2000)) = 0;
train_label=double(label(Indices~=j,1));  %%训练集标签

test_data2000=data2000(:,Indices==j)';     %%测试集数据
test_data2000(isnan(test_data2000)) = 0;
test_label=double(label(Indices==j,1));  %%测试集标签

nb_svm=fitcsvm(train_data2000,train_label);    %%SVM的训练模型
[res_svm, score] = predict(nb_svm,test_data2000);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分

[~,~,~,auc]=perfcurve(test_label,score(:,2),'1'); %计算分类器的AUC

auc_svm=auc;
cpre=classperf(test_label,res_svm);  %计算分类器的准确率等指标
Accuracy_svm=cpre.CorrectRate;  %%准确率
sensitivity_svm=cpre.Sensitivity;  %%灵敏度
specificity_svm=cpre.Specificity;   %%特异度
result=[auc_svm Accuracy_svm sensitivity_svm specificity_svm]./1;

title={'AUC','Accuracy','Seneitivity','Specificity'};
res=[title;num2cell(result)];
ressavename = strcat(dqpath,fname);
ressavename = strcat(ressavename,'_2000.mat');
save(ressavename,'res');
end