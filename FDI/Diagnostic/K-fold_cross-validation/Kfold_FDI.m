function Kfold_FDI(j,c1,h1,Indices,dqpath,cfdnapath,X,Y,readlen,buchang,map_path,gp,fdr)%优化了合并速度，但是算特征很慢，可以判断peak存在

fname = num2str(j);%What fold, the whole folder, dqpath represents the folder for this fold
mkdir(fname);
jk_dirname = strcat(fname,'/jk_traindata/');
mkdir(jk_dirname);
hcc_dirname = strcat(fname,'/hcc_traindata/');
mkdir(hcc_dirname);
dqpath = strcat(dqpath,fname);
dqpath = strcat(dqpath,'/');


x = 20; % Suppose you need to select 20 rows
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



dataname = [c1;h1];%



    
    for i =1:length(dataname)
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
    end
    for i =1:22
        ii = num2str(i);
        name = strcat('chr',ii);
        read = zeros(1,2);
        save([hcc_dirname, name],'read');
    end
    parfor i7 = 1:22
        chrname = strcat('chr',num2str(i7));
        chrname = strcat(chrname,'.mat');
        l3 = strcat(jk_dirname,chrname);
        s = load(l3);
        read_0 = s.read;
        for i1 = 1:x
            loadname = cell2mat(c2(i1));
            loadname = strcat(cfdnapath,loadname);
            loadname = strcat(loadname,'/data_n/');
            
            loadname1 = strcat(loadname,chrname );
            
            s1 = load(loadname1);
            read_0 = [read_0;s1.read];
        end
        read_0(~any(read_0,2),:)=[];
        read_0 = sortrows(read_0);
        read = read_0;
        savename = jk_dirname;
        r_parsave(savename,chrname,read,s1.region_len,s1.dark_flag);
        
    end
    parfor i7 = 1:22%
        chrname = strcat('chr',num2str(i7));
        chrname = strcat(chrname,'.mat');
        l3 = strcat(hcc_dirname,chrname);
        s = load(l3);
        read_0 = s.read;
        for i1 = 1:x
            loadname = cell2mat(h2(i1));
            loadname = strcat(cfdnapath,loadname);
            loadname = strcat(loadname,'/data_n/');
            
            loadname1 = strcat(loadname,chrname );
            
            s2 = load(loadname1);
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
    parfor i5 =1:2
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
            if peak_a(i+1,2)<peak_a(i,3) && peak_a(i+1,1)==peak_a(i,1)
                peak_a(i+1,2)=peak_a(i,2);
                peak_a(i+1,3) = max(peak_a(i+1,3),peak_a(i,3));
                peak_a(i+1,4)=max(peak_a(i+1,4),peak_a(i,4));
                peak_a(i+1,7)=max(peak_a(i+1,7),peak_a(i,7));
                peak_a(i+1,5)=max(peak_a(i+1,5),peak_a(i,5));
                peak_a(i+1,6)=max(peak_a(i+1,6),peak_a(i,6));
                peak_a(i+1,8)=peak_a(i+1,3)-peak_a(i+1,2);
                peak_a(i,:)=0;
            end
        end
    end
    peak_a(~any(peak_a,2),:)=[];
    peak_a = sortrows(peak_a,1);
    save([dqpath,'peak_all'],'peak_a');
    

parfor ij = 1:length(dataname)
    
    i55 = cell2mat(dataname(ij));
    
    feature2000 = getfeature2000(i55,dqpath,cfdnapath,X,Y);
    ln2 = strcat(dqpath,i55);
    ln2 = strcat(ln2,'/');
    parsave(feature2000(:,1),ln2,'feature');
    parsave(feature2000(:,2),ln2,'feature2000');
end

label  = ismember(dataname,h1);
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
%
train_data=data(:,Indices~=j)';     
train_data(isnan(train_data)) = 0;
train_label=double(label(Indices~=j,1)); 

test_data=data(:,Indices==j)';    
test_data(isnan(test_data)) = 0;
test_label=double(label(Indices==j,1));

nb_svm=fitcsvm(train_data,train_label);   
[res_svm, score] = predict(nb_svm,test_data);       

[~,~,~,auc]=perfcurve(test_label,score(:,2),'1');

auc_svm=auc;
cpre=classperf(test_label,res_svm);  
Accuracy_svm=cpre.CorrectRate;  
sensitivity_svm=cpre.Sensitivity; 
specificity_svm=cpre.Specificity; 
result=[auc_svm Accuracy_svm sensitivity_svm specificity_svm]./1;

title={'AUC','Accuracy','Seneitivity','Specificity'};
res=[title;num2cell(result)];
ressavename = strcat(dqpath,fname);
ressavename = strcat(ressavename,'.mat');
save(ressavename,'res');
train_data2000=data2000(:,Indices~=j)';     %%
train_data2000(isnan(train_data2000)) = 0;
train_label=double(label(Indices~=j,1));  %%

test_data2000=data2000(:,Indices==j)';     %%
test_data2000(isnan(test_data2000)) = 0;
test_label=double(label(Indices==j,1));  %%

nb_svm=fitcsvm(train_data2000,train_label);    %%
[res_svm, score] = predict(nb_svm,test_data2000);          %%
[~,~,~,auc]=perfcurve(test_label,score(:,2),'1'); %

auc_svm=auc;
cpre=classperf(test_label,res_svm);  %
Accuracy_svm=cpre.CorrectRate;  %
sensitivity_svm=cpre.Sensitivity;  %%
specificity_svm=cpre.Specificity;   %%
result=[auc_svm Accuracy_svm sensitivity_svm specificity_svm]./1;

title={'AUC','Accuracy','Seneitivity','Specificity'};
res=[title;num2cell(result)];
ressavename = strcat(dqpath,fname);
ressavename = strcat(ressavename,'_2000.mat');
save(ressavename,'res');
end