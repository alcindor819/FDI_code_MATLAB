
function classify_mutil(sample_info_path,name_col,fold,mat_folder)

    cfdnapath = mat_folder;
        [~,~,r] = xlsread(sample_info_path);
   
    
    h1 = r(ismember(r(:,name_col),'Breast_Cancer'),:);
    h2 = r(ismember(r(:,name_col),'Cholangiocarcinoma'),:);
    h3 = r(ismember(r(:,name_col),'Colorectal_Cancer'),:);
    h4 = r(ismember(r(:,name_col),'Gastric_cancer'),:);
    h5 = r(ismember(r(:,name_col),'Lung_Cancer'),:);
    h6 = r(ismember(r(:,name_col),'Ovarian_Cancer'),:);
    h7 = r(ismember(r(:,name_col),'Pancreatic_Cancer'),:);

    h1 = h1(:,1);
    h2 = h2(:,1);
    h3 = h3(:,1);
    h4 = h4(:,1);       
    h5 = h5(:,1);
    h6 = h6(:,1);
    h7 = h7(:,1);
    dataname = [h1;h2;h3;h4;h5;h6;h7];
    dqpath=pwd;%当前文件夹路径
    dqpath = strcat(dqpath,'/');
    rng('default');
    N = length(dataname);%RNG DEFULT
    n1 = zeros(N,1);
    for ii = 1:size(dataname,1)
        if ismember(dataname(ii),h1)
            n1(ii) = 1;
        end
        if ismember(dataname(ii),h2)
            n1(ii) = 2;
        end
        if ismember(dataname(ii),h3)
            n1(ii) = 3;
        end
        if ismember(dataname(ii),h4)
            n1(ii) = 4; 
        end
        if ismember(dataname(ii),h5)
            n1(ii) = 5;
        end
        if ismember(dataname(ii),h6)
            n1(ii) = 6;
        end
        if ismember(dataname(ii),h7)
            n1(ii) = 7;
        end
    end
    Indices = crossvalind('Kfold', n1, fold);   %%交叉验证，将样本分成3份

for j =1:10
    
    kfold_mutil(j,h1,h2,h3,h4,h5,h6,h7,Indices,dqpath,cfdnapath,n1);
end

dqpath = strcat(pwd,'/');
result = zeros(7,8);

for i =1:10
    name = num2str(i);
    d1 = strcat(dqpath,name);
    d2 = strcat(d1,'/');
    d2 = strcat(d2,num2str(i));
    d3 = strcat(d2,'.mat');
    d4 = strcat(d2,'_2000.mat');
    load(d3);
    
    for i3  = 1:7
        
        for i4 = 1:size(test_label,1)
        
                if  test_label(i4)==i3
                    result(i3,3) = result(i3,3)+1;
                    
                    if (test_label(i4) == res(i4))
                        result(i3,1) = result(i3,1)+1;
                    end
                
                end
        end
    end
    load(d4);
    for i3  = 1:7
        
        for i4 = 1:size(test_label,1)
        
                if  test_label(i4)==i3
                    result(i3,4) = result(i3,4)+1;
                    
                    if (test_label(i4) == res(i4))
                        result(i3,2) = result(i3,2)+1;
                    end
                
                end
        end
    end
    for i3  = 1:7
        
        for i4 = 1:size(test_label,1)
        
                if  secondLabels(i4)==i3
                    result(i3,7) = result(i3,7)+1;
                    
                    if (secondLabels(i4) == res(i4))
                        result(i3,5) = result(i3,5)+1;
                    end
                
                end
        end
    end
    for i3  = 1:7
        
        for i4 = 1:size(test_label,1)
        
                if  secondLabels2(i4)==i3
                    result(i3,8) = result(i3,8)+1;
                    
                    if (secondLabels2(i4) == res(i4))
                        result(i3,6) = result(i3,6)+1;
                    end
                
                end
        end
    end
end
    
    

save([dqpath,'mulclass_res'],'result');

end


