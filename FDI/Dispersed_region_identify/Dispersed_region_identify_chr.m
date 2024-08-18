% data_name = 'D:\wyzwork\Basic_info\IH02_data\';
% distribution = 'Beta';
% P = 0.01;
% FDR = 0.05;
% res_path = 'D:\wyzwork\精简版代码\FDI\Dispersed_region_identify\call\';
% res_name = 'IH02_test';
% X = 0.5;
% Y = 10;
% chr_n = 1;
function Dispersed_region_identify_chr(data_name,distribution,P,FDR,res_path,X,Y,res_name,chr_n)
result_name = res_name;
m2 = X;
peak_a=zeros(1,9);
Interval_length = 200;
footsteps = 20;
connection_parameter=200;
thresholdp=P;
thresholdfdr=FDR;
rb = Y;
file_path = data_name;
fd0 = 100;

dqpath = pwd;
positions = strfind(dqpath, '\');
lastPosition = positions(end);
lastpath = dqpath(1:lastPosition);
bipath = strcat(lastpath,'Basic_info\');

izong = chr_n;
izongchar=num2str(izong);

te=strcat('\data_n\chr',izongchar);
teee=strcat(file_path,te);
te1=strcat(teee,'.mat');
load (te1); 



file_path2=strcat(bipath,'mappability\');
ma1=strcat('chr',izongchar);
ma2=strcat(file_path2,ma1);
ma2 = strcat(ma2,'_m');
ma3=strcat(ma2,'.mat');
load (ma3);

read(:,3) = read(:,1)+read(:,2);
coverage = zeros(region_len,1);
for ii=1:length(read)
    for jj=read(ii,1):read(ii,1)+read(ii,2)
        if jj>region_len-1
            jj = region_len-1;
        end
        coverage(jj)=coverage(jj)+1;
    end
end
midu = zeros(region_len,2);
for i2 = 1:length(read)
    qist = read(i2,1);
    if qist<rb+1
        qist = rb+1;
    end
    if qist+rb>region_len-1
        qist = region_len-rb-1;
    end
    zhnd = read(i2,3);
    if zhnd<rb+1
        zhnd = rb+1;
    end
    if zhnd+rb>region_len-1
        zhnd = region_len-rb-1;
    end
    midu(qist,1) = midu(qist,1)+1;
    midu(zhnd,1) = midu(zhnd,1)+1;
    
end
for i2 = 1:length(read)
    qist = read(i2,1);
    if qist<rb+1
        qist = rb+1;
    end
    if qist+rb>region_len-1
        qist = region_len-rb-1;
    end
    zhnd = read(i2,3);
    if zhnd<rb+1
        zhnd = rb+1;
    end
    if zhnd+rb>region_len-1
        zhnd = region_len-rb-1;
    end
    zs = sum(midu((qist-rb:qist+rb),1))-1;
    midu(qist,2) = midu(qist,2)+m2^zs;
    zs = sum(midu((zhnd-rb:zhnd+rb),1))-1;
    midu(zhnd,2) = midu(zhnd,2)+m2^zs;
end



m=1;
countm = 0;
count_0 = 0;
count_reads = floor(length(coverage)\footsteps)+1;
dispersed_regions= zeros(count_reads,7);
for i=1:count_reads
    for k=m:m+Interval_length-1
        if k>(region_len-1)
            k=(region_len-1);
        end
        count_0 = count_0+coverage(k);
        if dark_flag(k) == 1
            countm=countm+1;
        end
    end
    mx1 = m;
    mx2 = m+Interval_length-1;
    if mx2>(region_len-1)
        mx2=(region_len-1);
    end
    countx = mean(map_b(mx1:mx2));
    
    iii=1;
    c_c = zeros(Interval_length,1);
    
    if countm == Interval_length && count_0 ~= 0 && countx>0.9
        
        fm200 = sum(midu(mx1:mx2,1));
        if fm200 ==0
            fm200 = 1;
        end
        edi = sum(midu(mx1:mx2,2))\fm200;
        
        
        
        for kk = m:m+Interval_length-1
            if kk >(length(coverage)-1)
                kk=(length(coverage)-1);
            end
            c_c(iii) = coverage(kk);
            iii = iii+1;
        end
        
        
        iii = 1;
        std_cov = std(c_c);
        
        
        
        dispersed_regions(i,6) = edi;
        
        dispersed_regions(i,3) = std_cov;
        dispersed_regions(i,1) = m;
        dispersed_regions(i,2) = m+Interval_length-1;
    end
    m=m+footsteps;
    countm=0;
    count_0=0;
    if m>length(coverage)
        m=length(coverage);
    end
end

for o6 =1:length(dispersed_regions)
    dispersed_regions(o6,3) = dispersed_regions(o6,3)*(dispersed_regions(o6,6));
end
dispersed_regions(:,4) = 0;
dispersed_regions(:,5) = 0;
dispersed_regions(:,6) = 0;
dispersed_regions(~any(dispersed_regions,2),:)=[];
data_true=dispersed_regions(:,3);
data_true = data_true\max(data_true);
pd = fitdist(data_true,distribution);

 %   data_expect=random(distribution,pd.a,pd.b,length(data_true),1);
pvalue = zeros(length(dispersed_regions),1);
%     p_expect=1-cdf(pd,data_expect);
% p_true=1-cdf(pd,data_true);
% 
% qqplot(-log10(p_expect),-log10(p_true));
% title('QQ plot');
% xlabel('Expected (-log10(P-value))');
% ylabel('Observed (-log10(P-value))');






for i=1:length(dispersed_regions)
    pvalue(i) = 1-cdf(pd,data_true(i));
end
for i=1:length(dispersed_regions)
    dispersed_regions(i,4) = pvalue(i);
end
FDR = mafdr(pvalue);
for i=1:length(dispersed_regions)
    dispersed_regions(i,5) = FDR(i);
end

fivekb = zeros(250,1);
count = 1;
ccc=0;
for i = 1:250
    fivekb(i) = 1;
end
for i=125:size(dispersed_regions,1)-125
    jjjr=125;
    jjjl=1;
    if (dispersed_regions(i,5)<thresholdfdr)&&(dispersed_regions(i,4)<thresholdp)
        ccc=ccc+1;
        for j=i:i+125
            fivekb(jjjr)=dispersed_regions(j,3);
            jjjr=jjjr+1;
        end
        for j= i-124:i-1
            fivekb(jjjl)=dispersed_regions(j,3);
            jjjl=jjjl+1;
        end
        pd1=fitdist(fivekb,'normal');%局部统计信息
        dispersed_regions(i,6) = (1-cdf(pd1,fivekb(125)));
        if dispersed_regions(i,6)<thresholdp
            dispersed_regions(i,7)=1;%作为一个暂时的标签，等会判断是否满足显著
        end
    end
end
chr=izong;
nr=sum(dispersed_regions(:,7));
notable_region = zeros(nr,7);
j=1;
for i=1:length(dispersed_regions)
    if dispersed_regions(i,7)==1
        notable_region(j,:)=dispersed_regions(i,:);
        j=j+1;
    end
end
notable_region(~any(notable_region,2),:)=[];
notable_region(:,9) = notable_region(:,2)-notable_region(:,1)+1;
co = 0;
for i=1:length(notable_region)
    
    notable_region(i,8) = co;
    co = 0;
end
notable_region(~any(notable_region,2),:)=[];

data_true6 = notable_region(:,8);
fd2 = prctile(data_true6, fd0);
for i=1:length(data_true6)
    
    if notable_region(i,8)>fd2
        notable_region(i,:) =0;
    end
end
notable_region(~any(notable_region,2),:)=[];
mm = size(notable_region,1);
if mm>2
    for i=1:size(notable_region,1)-1
        if notable_region(i+1,1)<notable_region(i,2)%连接一下互相重叠的
            notable_region(i+1,1)=notable_region(i,1);
            notable_region(i+1,3)=max(notable_region(i+1,3),notable_region(i,3));
            notable_region(i+1,4)=max(notable_region(i+1,4),notable_region(i,4));
            notable_region(i+1,5)=max(notable_region(i+1,5),notable_region(i,5));
            notable_region(i+1,6)=max(notable_region(i+1,6),notable_region(i,6));
            notable_region(i,:)=0;
        end
    end
end
notable_region(~any(notable_region,2),:)=[];

for i=1:size(notable_region,1)-1
    if (notable_region(i+1,1)-notable_region(i,2))<connection_parameter%连接一下起点很近的
        notable_region(i+1,1)=notable_region(i,1);
        notable_region(i+1,3)=max(notable_region(i+1,3),notable_region(i,3));
        notable_region(i+1,4)=max(notable_region(i+1,4),notable_region(i,4));
        notable_region(i+1,5)=max(notable_region(i+1,5),notable_region(i,5));
        notable_region(i+1,6)=max(notable_region(i+1,6),notable_region(i,6));
        notable_region(i,:)=0;
    end
end

notable_region(~any(notable_region,2),:)=[];

mmm = size(notable_region,1);

chh = zeros(mmm,1);
for i=1:mmm
    chh(i)=chr;
end
n_r = [notable_region chh];
teee=strcat(te,'.mat');
fnr = zeros(size(n_r,1),8);
fnr(:,1) = n_r(:,10);
fnr(:,2) = n_r(:,1);
fnr(:,3) = n_r(:,2);
fnr(:,4) = n_r(:,3);
fnr(:,5) = n_r(:,4);
fnr(:,6) = n_r(:,5);
fnr(:,7) = n_r(:,6);
fnr(:,8) = fnr(:,3)-fnr(:,2)+1;
fnr(:,9) = n_r(:,8);
peak_a=[peak_a;fnr];

peak_a(~any(peak_a,2),:)=[];
peak_a = sortrows(peak_a);
peak_a(:,9) = [];
lujing1 = res_path;
%mkdir (lujing1);
lujing2 = strcat(lujing1,'\');
lujing2 = strcat(lujing2,result_name);
luj = '\result_n';

lujing3 = strcat(lujing2,luj);
mkdir (lujing3);
lujing4 = strcat(lujing3,'\');
save([lujing4, strcat('peak_',num2str(chr_n))],'peak_a');
flag = strcat('chr_',num2str(izong));
flag = strcat(flag,' has been finished');
display(flag);


end