function feature2 = getfeature2000(name,dqpath,cfdnapath,X,Y)%修正端点离散度
i55 = name;
m2 = X;
rb = Y;
feature = zeros(1,2);
l3 = strcat(dqpath,'peak_all.mat');
load(l3);
for i6  = 1:22
    i66 = strcat('chr',num2str(i6));
    i66 = strcat('/data_n/',i66);
    i66 = strcat(i55,i66);
    ln = strcat(cfdnapath,i66);
    ln = strcat(ln,'.mat');
    load(ln);
    peak_a1 = peak_a(ismember(peak_a(:,1),i6),:);
    if isempty(peak_a1)
        continue;
    end
    coverage = zeros(region_len,1);
    read(:,3) = read(:,1)+read(:,2);
    for ii=1:length(read)
        for jj=read(ii,1):read(ii,1)+read(ii,2)
            if jj>region_len-1
                jj = region_len-1;
            end
            coverage(jj)=coverage(jj)+1;
        end
    end%计算当前染色体覆盖度
    midu = zeros(region_len,2);
    for i7 = 1:length(read)
        qist = ceil(read(i7,1));
        if qist<rb+1
            qist = rb+1;
        end
        if qist+rb>region_len-1
            qist = region_len-rb-1;
        end
        zhnd = ceil(read(i7,3));
        if zhnd<rb+1
            zhnd = rb+1;
        end
        if zhnd+rb>region_len-1
            zhnd = region_len-rb-1;
        end
        midu(qist,1) = midu(qist,1)+1;
        midu(zhnd,1) = midu(zhnd,1)+1;

    end
    for i7 = 1:length(read)
        qist = read(i7,1);
        if qist<rb+1
            qist = rb+1;
        end
        if qist+rb>region_len-1
            qist = region_len-rb-1;
        end
        zhnd = read(i7,3);
        if zhnd<rb+1
            zhnd = rb+1;
        end
        if zhnd+rb>region_len-1
            zhnd = region_len-rb-1;
        end
        bbb1 = ceil(qist-rb);
        bbb2 = ceil(qist+rb);
        bbb3 = ceil(zhnd-rb);
        bbb4 = ceil(zhnd+rb);
        zs = sum(midu((bbb1:bbb2),1))-1;
        
        midu(qist,2) = midu(qist,2)+m2^zs;
        zs = sum(midu((bbb3:bbb4),1))-1;
        midu(zhnd,2) = midu(zhnd,2)+m2^zs;
    end%计算当前染色体的端点离散度
    feature1 = zeros(size(peak_a1,1),2);
    for i8 = 1:size(feature1,1)
        mx0 = ceil(0.5*(peak_a1(i8,2)+peak_a1(i8,3)));
        mx1 = mx0-1000;
        mx2 = mx0+1000;
        mx11 = peak_a1(i8,2);
        mx22 = peak_a1(i8,3);
        if mx1<1
            mx1 = 1;
        end
        if mx2>region_len
            mx2 = region_len;
        end
        c2000 = std(coverage(mx1:mx2));
        c = std(coverage(mx11:mx22));
        fm2000 = sum(midu(mx1:mx2,1));
        if fm2000 ==0
            fm2000 = 1;
        end
        wyzp2000 = sum(midu(mx1:mx2,2))/fm2000;
        p2000 = wyzp2000;
        
        fm = sum(midu(mx11:mx22,1));
        if fm ==0
            fm = 1;
        end
        wyzp = sum(midu(mx11:mx22,2))/fm;
        p = wyzp;
        feature1(i8,1) = c*p;
        feature1(i8,2) = c2000*p2000;
    end%计算当前染色体上的特征
    feature = [feature;feature1];
    clear feature1;
end
feature(1,:) = [];
feature2 = feature;
end