function kfold_mutil(j,h1,h2,h3,h4,h5,h6,h7,Indices,dqpath,cfdnapath,n1)%下采样的多分类,保存预测结果

positions = strfind(dqpath, '/');
lastPosition = positions(end-2);
lastpath = dqpath(1:lastPosition);
bipath = strcat(lastpath,'Basic_info/');
map_path = strcat(bipath,'mappability/');

fname = num2str(j);
mkdir(fname);

h1_dirname = strcat(fname,'/h1_traindata/');
mkdir(h1_dirname);
h2_dirname = strcat(fname,'/h2_traindata/');
mkdir(h2_dirname);
h3_dirname = strcat(fname,'/h3_traindata/');
mkdir(h3_dirname);
h4_dirname = strcat(fname,'/h4_traindata/');
mkdir(h4_dirname);
h5_dirname = strcat(fname,'/h5_traindata/');
mkdir(h5_dirname);
h6_dirname = strcat(fname,'/h6_traindata/');
mkdir(h6_dirname);
h7_dirname = strcat(fname,'/h7_traindata/');
mkdir(h7_dirname);


dqpath = strcat(dqpath,fname);
dqpath = strcat(dqpath,'/');
dataname = [h1;h2;h3;h4;h5;h6;h7];
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
    
    save([h1_dirname, name],'read');
    save([h2_dirname, name],'read');
    save([h3_dirname, name],'read');
    save([h4_dirname, name],'read');
    save([h5_dirname, name],'read');
    save([h6_dirname, name],'read');
    save([h7_dirname, name],'read');
end%生成初始chr

train_name = dataname(Indices~=j);%综合训练集名字
test_name = dataname(Indices==j);%综合测试集名字
h1tname = train_name(ismember(train_name,h1));%HCC训练集名字
h2tname = train_name(ismember(train_name,h2));%HCC训练集名字
h3tname = train_name(ismember(train_name,h3));%HCC训练集名字
h4tname = train_name(ismember(train_name,h4));%HCC训练集名字
h5tname = train_name(ismember(train_name,h5));%HCC训练集名字
h6tname = train_name(ismember(train_name,h6));%HCC训练集名字
h7tname = train_name(ismember(train_name,h7));%HCC训练集名字
sumtrain = 20;

parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h1_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:sumtrain%整合JK训练集的CFDNA
        loadname = cell2mat(h1tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s1 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s1.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h1_dirname;
    r_parsave(savename,chrname,read,s1.region_len,s1.dark_flag);
    
end
parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h2_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:sumtrain%整合JK训练集的CFDNA
        loadname = cell2mat(h2tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s2 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s2.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h2_dirname;
    r_parsave(savename,chrname,read,s2.region_len,s2.dark_flag);
    
end
parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h3_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:sumtrain%整合JK训练集的CFDNA
        loadname = cell2mat(h3tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s3 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s3.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h3_dirname;
    r_parsave(savename,chrname,read,s3.region_len,s3.dark_flag);
    
end
parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h3_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:sumtrain%整合JK训练集的CFDNA
        loadname = cell2mat(h4tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s4 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s4.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h4_dirname;
    r_parsave(savename,chrname,read,s4.region_len,s4.dark_flag);
    
end    
parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h3_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:8%整合JK训练集的CFDNA,肺癌的少
        loadname = cell2mat(h5tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s5 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s5.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h5_dirname;
    r_parsave(savename,chrname,read,s5.region_len,s5.dark_flag);
    
end
parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h3_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:sumtrain%整合JK训练集的CFDNA
        loadname = cell2mat(h6tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s6 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s6.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h6_dirname;
    r_parsave(savename,chrname,read,s6.region_len,s6.dark_flag);
    
end
parfor i7 = 1:22%把新的CFDNA加进去
    chrname = strcat('chr',num2str(i7));
    chrname = strcat(chrname,'.mat');
    l3 = strcat(h3_dirname,chrname);
    s = load(l3);
    read_0 = s.read;
    for i1 = 1:sumtrain%整合JK训练集的CFDNA
        loadname = cell2mat(h7tname(i1));
        loadname = strcat(cfdnapath,loadname);
        loadname = strcat(loadname,'/data_n/');
        
        loadname1 = strcat(loadname,chrname );
        
        s7 = load(loadname1);%读取将要合并的数据
        read_0 = [read_0;s7.read];
    end
    read_0(~any(read_0,2),:)=[];
    read_0 = sortrows(read_0);
    read = read_0;
    savename = h7_dirname;
    r_parsave(savename,chrname,read,s7.region_len,s7.dark_flag);
    
end
parfor i5 =2:8
    if i5==1
        jkjkjk = 1;
    elseif i5 ==2
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h1peak',0.5,10,200,20,map_path);
    elseif i5 ==3
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h2peak',0.5,10,200,20,map_path);
    elseif i5 ==4
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h3peak',0.5,10,200,20,map_path);
    elseif i5 ==5
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h4peak',0.5,10,200,20,map_path);
    elseif i5 ==6
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h5peak',0.5,10,200,20,map_path);
    elseif i5 ==7
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h6peak',0.5,10,200,20,map_path);
    elseif i5 ==8
        FDI_call(0.001,0.05,h1_dirname,h1_dirname,'h7peak',0.5,10,200,20,map_path);
        
        
    end
end



hebingh1peak = strcat(h1_dirname,'h1peak/result_n/peak_all.mat');
load(hebingh1peak);
peak_a1 = peak_a;


peak_a1 = peak_a;
hebingh2peak = strcat(h2_dirname,'h2peak/result_n/peak_all.mat');
load(hebingh2peak);
peak_a = [peak_a1;peak_a];

peak_a1 = peak_a;
hebingh3peak = strcat(h3_dirname,'h3peak/result_n/peak_all.mat');
load(hebingh3peak);
peak_a = [peak_a1;peak_a];

peak_a1 = peak_a;
hebingh4peak = strcat(h4_dirname,'h4peak/result_n/peak_all.mat');
load(hebingh4peak);
peak_a = [peak_a1;peak_a];

peak_a1 = peak_a;
hebingh5peak = strcat(h5_dirname,'h5peak/result_n/peak_all.mat');
load(hebingh5peak);
peak_a = [peak_a1;peak_a];

peak_a1 = peak_a;
hebingh6peak = strcat(h6_dirname,'h6peak/result_n/peak_all.mat');
load(hebingh6peak);
peak_a = [peak_a1;peak_a];

peak_a1 = peak_a;
hebingh7peak = strcat(h7_dirname,'h7peak/result_n/peak_all.mat');
load(hebingh7peak);
peak_a = [peak_a1;peak_a];

peak_a = sortrows(peak_a,[1 2]);



mm = size(peak_a,1);
if mm>2
    for i=1:size(peak_a,1)-1
        if peak_a(i+1,2)<peak_a(i,3) && peak_a(i+1,1)==peak_a(i,1)%连接一下互相重叠的
            peak_a(i+1,2)=peak_a(i,2);
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
parfor i5 = 1:length(dataname)
    
    i55 = cell2mat(dataname(i5));%具体人的文件夹名字
    
    feature2000 = getfeature2000(i55,dqpath,cfdnapath);
    ln2 = strcat(dqpath,i55);
    ln2 = strcat(ln2,'/');
    parsave(feature2000(:,1),ln2,'feature');
    parsave(feature2000(:,2),ln2,'feature2000');
end

label  = n1;
for i9 =1:length(dataname)
    ln = strcat(dqpath,cell2mat(dataname(i9)));
    ln = strcat(ln,'/feature.mat');
    load(ln);
    feature_all(:,i9) = feature3;%一列是一个样本，一行对应一个区间，
end

for i6 = 1:size(feature_all,1)
    
    
    C1 = zeros(length(c1),1);
    i666 = 1;
    for i66 = 1:length(c1)
        C1(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g0=repmat({'C1'},size(C1,1));
    
    H1 = zeros(length(h1),1);
    i666 = 1;
    for i66 = 0+1:0+length(h1)
        H1(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g1=repmat({'H1'},size(H1,1));
    
    H2 = zeros(length(h2),1);
    i666 = 1;
    for i66 = 0+length(h1)+1:0+length(h1)+length(h2)
        H2(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g2=repmat({'H2'},size(H2,1));
    
    H3 = zeros(length(h3),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+1:0+length(h1)+length(h2)+length(h3)
        H3(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g3=repmat({'H3'},size(H3,1));
    
    H4 = zeros(length(h4),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+1:0+length(h1)+length(h2)+length(h3)+length(h4)
        H4(i666) = feature_dall(i6,i66);
        i666 = i666+1;
    end
    g4=repmat({'H4'},size(H4,1));
    
    H5 = zeros(length(h5),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+length(h4)+1:0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)
        H5(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g5=repmat({'H5'},size(H5));
    
    H6 = zeros(length(h6),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+1:0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+length(h6)
        H6(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g6=repmat({'H6'},size(H6,1));
    
    H7 = zeros(length(h7),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+length(h6)+1:0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+length(h6)+length(h7)
        H7(i666) = feature_all(i6,i66);
        i666 = i666+1;
    end
    g7=repmat({'H7'},size(H7,1));
    life=[H1;H2;H3;H4;H5;H6;H7];
    %g0 = g0(:,1);
    
    
    g1 = g1(:,1);
    g2 = g2(:,1);
    g3 = g3(:,1);
    g4 = g4(:,1);
    g5 = g5(:,1);
    g6 = g6(:,1);
    g7 = g7(:,1);
    group=[g1;g2;g3;g4;g5;g6;g7];


    p2=anova1(life,group,'off');
    if p2>0.05
        feature_all(i6,:) = 0;
    end
    %如果一行没通过检验，直接删除
end

for i9 =1:length(dataname)
    ln = strcat(dqpath,cell2mat(dataname(i9)));
    ln = strcat(ln,'/feature2000.mat');
    load(ln);
    feature_all2000(:,i9) = feature3;
end
for i6 = 1:size(feature_all,1)
    
    
    H1 = zeros(length(h1),1);
    i666 = 1;
    for i66 = 0+1:0+length(h1)
        H1(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g1=repmat({'H1'},size(H1,1));
    
    H2 = zeros(length(h2),1);
    i666 = 1;
    for i66 = 0+length(h1)+1:0+length(h1)+length(h2)
        H2(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g2=repmat({'H2'},size(H2,1));
    
    H3 = zeros(length(h3),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+1:0+length(h1)+length(h2)+length(h3)
        H3(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g3=repmat({'H3'},size(H3,1));
    
    H4 = zeros(length(h4),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+1:0+length(h1)+length(h2)+length(h3)+length(h4)
        H4(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g4=repmat({'H4'},size(H4,1));
    
    H5 = zeros(length(h5),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+length(h4)+1:0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)
        H5(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g5=repmat({'H5'},size(H5));
    
    H6 = zeros(length(h6),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+1:0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+length(h6)
        H6(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g6=repmat({'H6'},size(H6,1));
    
    H7 = zeros(length(h7),1);
    i666 = 1;
    for i66 = 0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+length(h6)+1:0+length(h1)+length(h2)+length(h3)+length(h4)+length(h5)+length(h6)+length(h7)
        H7(i666) = feature_all2000(i6,i66);
        i666 = i666+1;
    end
    g7=repmat({'H7'},size(H7,1));
    life=[H1;H2;H3;H4;H5;H6;H7];
    
    g1 = g1(:,1);
    g2 = g2(:,1);
    g3 = g3(:,1);
    g4 = g4(:,1);
    g5 = g5(:,1);
    g6 = g6(:,1);
    g7 = g7(:,1);
    group=[g1;g2;g3;g4;g5;g6;g7];
    

    p2=anova1(life,group,'off');
    if p2>0.05
        feature_all2000(i6,:) = 0;
    end
    %如果一行没通过检验，直接删除
end




data = feature_all;
data2000 = feature_all2000;
data(isnan(data)) = 0;
data(all(data==0,2),:)=[];
data2000(isnan(data2000)) = 0;
data2000(all(data2000==0,2),:)=[];
% for i33 = 1:size(data,2)
%     data(:,i33) = zscore(data(:,i33));
% end
% for i33 = 1:size(data2000,2)
%     data2000(:,i33) = zscore(data2000(:,i33));
% end
%开始照搬分类器原始peak
train_data=data(:,Indices~=j)';     %%训练集数据
train_data(isnan(train_data)) = 0;
train_label=double(label(Indices~=j,1));  %%训练集标签

test_data=data(:,Indices==j)';     %%测试集数据
test_data(isnan(test_data)) = 0;
test_label=double(label(Indices==j,1));  %%测试集标签
X= train_data;
Y = train_label;
Z = test_data;

at1 = zeros(1,length(train_data));
at1(1,length(train_data)) = 1;
attribute = at1;
ClassType = [1,2,3,4,5,6,7];
C=[1,1,1,1,1,1,1];
addpath('CSNN/');
[X1,Y1]=SmoteOverSampling(X',Y',ClassType,C,attribute,5,'nominal');

Mdl = fitcecoc(X1',Y1');
Z1 = predict(Mdl,Z);



[~, score] = predict(Mdl, Z);

% 对于每个样本，找到最高和次高的分数及其对应的标签
[~, maxScoreIdx] = max(score, [], 2);
score(sub2ind(size(score), (1:size(score,1))', maxScoreIdx)) = -inf;
[~, secondMaxScoreIdx] = max(score, [], 2);

% 获取标签
labels = Mdl.ClassNames;
firstLabels = labels(maxScoreIdx);
secondLabels = labels(secondMaxScoreIdx);



s3 = 0;
for i3  = 1:size(test_label,1)
    if test_label(i3) == Z1(i3)
        s3 = s3+1;
    end
end
acc = s3/size(test_label,1);
res=Z1;
ressavename = strcat(dqpath,fname);
ressavename = strcat(ressavename,'.mat');
save(ressavename,'res','test_label');
%开始照搬分类器   2000peak
train_data2000=data2000(:,Indices~=j)';     %%训练集数据
train_data2000(isnan(train_data2000)) = 0;
train_label=double(label(Indices~=j,1));  %%训练集标签

test_data2000=data2000(:,Indices==j)';     %%测试集数据
test_data2000(isnan(test_data2000)) = 0;
test_label=double(label(Indices==j,1));  %%测试集标签

X= train_data2000;
Y = train_label;
Z = test_data2000;

at1 = zeros(1,length(train_data2000));
at1(1,length(train_data2000)) = 1;
attribute = at1;
ClassType = [1,2,3,4,5,6,7];
C=[1,1,1,1,1,1,1];

[X1,Y1]=SmoteOverSampling(X',Y',ClassType,C,attribute,5,'nominal');

Mdl = fitcecoc(X1',Y1');


Z1 = predict(Mdl,Z);

[~, score] = predict(Mdl, Z);

% 对于每个样本，找到最高和次高的分数及其对应的标签
[~, maxScoreIdx] = max(score, [], 2);
score(sub2ind(size(score), (1:size(score,1))', maxScoreIdx)) = -inf;
[~, secondMaxScoreIdx] = max(score, [], 2);

% 获取标签
labels = Mdl.ClassNames;
firstLabels2 = labels(maxScoreIdx);
secondLabels2 = labels(secondMaxScoreIdx);



s3 = 0;
for i3  = 1:size(test_label,1)
    if test_label(i3) == Z1(i3)
        s3 = s3+1;
    end
end
acc = s3/size(test_label,1);
res=Z1;
ressavename = strcat(dqpath,fname);
ressavename = strcat(ressavename,'_2000.mat');
save(ressavename,'res','test_label','firstLabels2','firstLabels','secondLabels','secondLabels2');
end