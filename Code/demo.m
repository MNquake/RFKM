clc,clear; %清除命令，清空工作区，关闭所有窗口
data = load('E:\Data\dataset\MSRA25.mat','X');
label = load('E:\Data\dataset\MSRA25.mat','Y');
data = struct2cell(data);
label = struct2cell(label);
data = cell2mat(data);
label = cell2mat(label);
data = double(data);
data = mapminmax(data,0,1);
data = data';
% label = label+1;
label = label';
k = 12;                        % 类别
K = 45;
cluster_name = 'Yale32';
num = 1;
i = 150;
alphat=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y'];
percent = 0;
excel = actxserver('Excel.Application');
workbook = excel.Workbooks.Open('E:\RFKM\result\1111.xlsx');
worksheet = workbook.Sheets.Item(1);
for r = 1.1 : 0.1 : 2 
for percent = 0 : 0.1 : 0.4
for K = 5 : 30 : 65
    begin = [alphat(1),num2str(num)];
    ends = [alphat(3),num2str(num)];
    ff = [begin,':',ends];
    parameter(1,1) = r;
    parameter(1,2) = percent;
    parameter(1,3) = K;
    range = worksheet.Range(ff);
    range.Value = parameter;
    workbook.Save();
%     xlswrite('G:\Data\RFKM\result.xlsx', parameter, 'Sheet1', ff);
    index = 1;
    tic
while 0<1
    if index == 26
        break;
    end
    [col, row]=size(data);
    U = rand(row, k);              %随机模糊矩阵
    row_sum= sum(U, 2);
    U = U./(row_sum*ones(1, k));   %约束条件：每一行累加为1
    F=U;
    
    S = ones(row, k);                   %噪声矩阵
    noise = fix(row * percent);         %假设20%为噪声
    
    zero_rows = randperm(row, noise);   %随机选择噪声点
    S(zero_rows,:) = 0;                 %选择为噪声点的位置置为0，表示为噪声，非噪声则为1
    
%     center = zeros(row,k);
    [F, obj_RFKM, iter, cid, S] = newRFKM(F, r, data, noise,K, S);
%     close all;
    [~,max_F] = max(F,[],2);
    [idx,~] = find(S(:,1)~=0);
    max_F = max_F(idx,:);
    label_ = label(idx,:);

%     [label_idx,~] = find(label_ == 0);
%     for t = label_idx
%         label_(t) = randi(k, 1);
%     end

    [Purity, ACC, ARI, NMI] = Evaluation(label_, max_F);
    
%     [data_idx, ~] = find(S(:,1) == 0);
%     [label_idx,~] = find(label == 0);
% 
%     C = intersect(data_idx,label_idx);
%     noiseACC = length(C) / length(data_idx);

    result(1,1)=iter;
    result(2,1)=obj_RFKM(iter);
    result(3,1)=Purity;
    result(4,1)=ARI;
    result(5,1)=ACC;
    result(6,1)=NMI;
%     result(7,1)=noiseACC;
    sss = alphat(index);
    begin = [sss,num2str(num + 1)];
    ends = [sss,num2str(num + 6)];
    fi = [begin,':',ends];
    range = worksheet.Range(fi);
    range.Value = result;
    workbook.Save();
    %xlswrite('G:\Data\RFKM\result.xlsx', result, 'Sheet2', fi);
    index = index + 1;
end
    num = num + 7;
    i = i - 1;
    time = (i * toc) / 3600;
    fprintf('此次迭代用时：%f s, 剩余迭代次数：%d 次，算法剩余时间：%f h\n',toc, i, time);
end
end
end
workbook.Close(false);
excel.Quit();

























