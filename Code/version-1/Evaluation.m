function [Purity, ACC, ARI, NMI] = Evaluation(label,result)

    NewLabel = BestMapping(label',result');      % 行向量
    
    Purity = Purity_cal(double(label'),NewLabel);

    ACC = ACC_cal(label',NewLabel);

    ARI = RandIndex(uint8(label'+1),uint8(result'+1));

    NMI = NMI_cal(label',result');


end

function [Purity] = Purity_cal(label,result)
    %label为人工标记簇
    %result为聚类算法计算结果
    N = length(result);%样本总数
    label_num = unique(label);
    k = unique(result);
    P_size = length(label_num); %人工标记的簇的个数
    C_size = length(k); %算法计算的簇的个数
    % Pid Cid 第i行非零数据代表的样本属于第i个簇
    Label_id = double(ones(P_size,1)*label==label_num'*ones(1,N));
    k_id = double(ones(C_size,1)*result==k'*ones(1,N));
    CP = k_id*Label_id'; %label和result的交集
    Pj = sum(CP,1);%label在result各个簇中的个数
    Ci = sum(CP,2);%result在lable各个簇中的个数
    
    precision = CP./(Ci*ones(1,P_size));
    recall = CP./(ones(C_size,1)*Pj);
    F = 2*precision.*recall./(precision+recall);
    FMeasure = sum((Pj./sum(Pj)).*max(F));
    Purity = sum(max(CP,[],2))/N;
end

function acc = ACC_cal(Label1,Label2)
    %Label1:真实标签 Label2:映射后的标签
    T= Label1==Label2;
    acc=sum(T)/length(Label2);
end


function [AR] = RandIndex(C1,C2)
    %C1为label
    %C2为聚类结果
    C = Contingency(C1,C2);
    n = sum(sum(C));
    nis = sum(sum(C,2).^2);
    njs = sum(sum(C,1).^2);
    t1 = nchoosek(n,2); %total number of pairs of entities
    t2 = sum(sum(C.^2)); %sum over rows & columnns of nij^2
    t3=.5*(nis+njs);
    nc = (n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));
    A = t1+t2-t3;
    if t1 == nc
        AR = 0;
    else
        AR=(A-nc)/(t1-nc);
    end
end

function Cont = Contingency(Mem1,Mem2)
    Cont = zeros(max(Mem1),max(Mem2));
    for i = 1:length(Mem1)
        Cont(Mem1(i),Mem2(i)) = Cont(Mem1(i),Mem2(i))+1;
    end
end

function NMI = NMI_cal(label, result)
    %NMI Normalized mutual information
    % A, B: 1*N;
    if length(label) ~= length(result)
        error('length( A ) must == length( B)');
    end
    N = length(label);
    A_id = unique(label);
    K_A = length(A_id);
    B_id = unique(result);
    K_B = length(B_id);
    % Mutual information
    A_occur = double (repmat( label, K_A, 1) == repmat( A_id', 1, N ));
    B_occur = double (repmat( result, K_B, 1) == repmat( B_id', 1, N ));
    AB_occur = A_occur * B_occur';
    P_A= sum(A_occur') / N;
    P_B =sum(B_occur') / N;
    P_AB = AB_occur / N;
    MImatrix = P_AB .* log(P_AB ./(P_A' * P_B)+eps);
    MI = sum(MImatrix(:));
    % Entropies
    H_A = -sum(P_A .* log(P_A + eps),2);
    H_B= -sum(P_B .* log(P_B + eps),2);
    %Normalized Mutual information
    NMI = MI / sqrt(H_A*H_B);
end

function [NewLabel] = BestMapping(La1,La2)
    %真实标签：La1 聚类结果标签：La2 映射后的标签：NewLabel

    Label1=unique(La1');
    L1=length(Label1);
    Label2=unique(La2');
    L2=length(Label2);

    %构建计算两种分类标签重复度的矩阵G
    G = zeros(max(L1,L2),max(L1,L2));
    for i=1:L1
        index1= La1==Label1(i,1);
        for j=1:L2
            index2= La2==Label2(j,1);
            G(i,j)=sum(index1.*index2);
        end
    end
    %利用匈牙利算法计算出映射重排后的矩阵
    [index]=munkres(-G);
    %将映射重排结果转换为一个存储有映射重排后标签顺序的行向量
    [temp]=MarkReplace(index);
    %生成映射重排后的标签NewLabel
    NewLabel=zeros(size(La2));
    for i=1:L2
        NewLabel(La2==Label2(i))=temp(i);
    end
end

function [assignment] = munkres(costMat)

    assignment = false(size(costMat));
     
    costMat(costMat~=costMat)=Inf;
    validMat = costMat<Inf;
    validCol = any(validMat);
    validRow = any(validMat,2);
     
    nRows = sum(validRow);
    nCols = sum(validCol);
    n = max(nRows,nCols);
    if ~n
        return
    end
         
    dMat = zeros(n);
    dMat(1:nRows,1:nCols) = costMat(validRow,validCol);
    dMat = bsxfun(@minus, dMat, min(dMat,[],2));
    zP = ~dMat;
    starZ = false(n);
    while any(zP(:))
        [r,c]=find(zP,1);
        starZ(r,c)=true;
        zP(r,:)=false;
        zP(:,c)=false;
    end
     
    while 1
        primeZ = false(n);
        coverColumn = any(starZ);
        if ~any(~coverColumn)
            break
        end
        coverRow = false(n,1);
        while 1
            zP(:) = false;
            zP(~coverRow,~coverColumn) = ~dMat(~coverRow,~coverColumn);
            Step = 6;
            while any(any(zP(~coverRow,~coverColumn)))
                [uZr,uZc] = find(zP,1);
                primeZ(uZr,uZc) = true;
                stz = starZ(uZr,:);
                if ~any(stz)
                    Step = 5;
                    break;
                end
                coverRow(uZr) = true;
                coverColumn(stz) = false;
                zP(uZr,:) = false;
                zP(~coverRow,stz) = ~dMat(~coverRow,stz);
            end
            if Step == 6
                M=dMat(~coverRow,~coverColumn);
                minval=min(min(M));
                if minval==inf
                    return
                end
                dMat(coverRow,coverColumn)=dMat(coverRow,coverColumn)+minval;
                dMat(~coverRow,~coverColumn)=M-minval;
            else
                break
            end
        end
        rowZ1 = starZ(:,uZc);
        starZ(uZr,uZc)=true;
        while any(rowZ1)
            starZ(rowZ1,uZc)=false;
            uZc = primeZ(rowZ1,:);
            uZr = rowZ1;
            rowZ1 = starZ(:,uZc);
            starZ(uZr,uZc)=true;
        end
    end
    assignment(validRow,validCol) = starZ(1:nRows,1:nCols);
end

%将存储标签顺序的空间矩阵转换为一个行向量
function [assignment] = MarkReplace(MarkMat)
    [rows,cols]=size(MarkMat);
    assignment=zeros(1,cols);
    for i=1:rows
        for j=1:cols
            if MarkMat(i,j)==1
                assignment(1,j)=i;
            end
        end
    end
end
