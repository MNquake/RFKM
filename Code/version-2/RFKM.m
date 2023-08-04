function [Z,e,obj_RFKM] = RFKM(d,k,h)
% This function aims to solve the Fuzzy k-Medoids problem:
% minimize Z = sum{i}sum{j}d{ij}*e{ij}^h
% subject to
% sum{j}e{ij} = 1 for all i
% sum{j}e{jj} = k
% e{ij} <= e{jj} for all i,j
% 0 <= e{ij} <= 1 for all i ~= j
% e{jj} binary for all j

% Inputs:
% d: n x n dissimilarity matrix
% k: number of clusters
% h: membership fuzziness factor

% Outputs:
% Z: solution cost
% e: k x n membership matrix

n = size(d,1);

% n个数据 k类
medoids = randperm(n,k);    %随机选择中心点
obj_RFKM = zeros(1, 100);   
uij = zeros(n,k);           %隶属矩阵
dij = d(:,medoids);         %所有点到中心点的距离             
if h == 1
    for i = 1:n
        [~,j] = min(dij(i,:));
        uij(i,j) = 1;
    end
else                                            %更新隶属矩阵
    for i = 1:n
        soma = 0;
        for j = 1:k
            soma = soma + 1/dij(i,j)^(1/(h-1));
        end
        for j = 1:k
            if dij(i,j) == 0
                uij(i,j) = 1;
            else
                uij(i,j) = 1/(soma*dij(i,j)^(1/(h-1)));
            end
        end
    end
end
Z = sum(sum(dij.*uij.^h));      %目标函数
obj_RFKM(1) = Z;
m_bin = zeros(n,1);             %记录中心点
for j = 1:k
    m_bin(medoids(j)) = 1;      %记录哪些点是中心点，哪些点不是
end

improved = true;
t = 1;
while improved
    improved = false;
    m = medoids;                %所有中心点
    m_bin2 = m_bin;
    for idOut = randperm(k)         %随机取中心点
        mOut = m(idOut);            %取出一个中心点
        for mIn = 1:n               %遍历点，用来更新中心点
            if m_bin2(mIn) == 1     
                continue;
            end
            m(idOut) = mIn;         % 先将这个点作为中心点
            m_bin2(mOut) = 0;       % 切换中心点
            m_bin2(mIn) = 1;
            dij = d(:,m);           % 更新其他点与中心点距离
            u2 = zeros(n,k);
            if h == 1               % 更新隶属矩阵
                for i = 1:n
                    [~,j] = min(dij(i,:));
                    u2(i,j) = 1;
                end
            else
                for i = 1:n
                    soma = 0;
                    for j = 1:k
                        soma = soma + 1/dij(i,j)^(1/(h-1));
                    end
                    for j = 1:k
                        if dij(i,j) == 0
                            u2(i,j) = 1;
                        else
                            u2(i,j) = 1/(soma*dij(i,j)^(1/(h-1)));
                        end
                    end
                end
            end
            Z2 = sum(sum(dij.*u2.^h));
            if Z2 < Z               % 此次更新使得函数目标值变小 记录更新
                t = t + 1;
                obj_RFKM(t) = Z2;
                improved = true;
                medoids = m;
                m_bin = m_bin2;
                Z = Z2;
                uij = u2;
            else
                m = medoids;
                m_bin2 = m_bin;
            end
        end
    end
end
e = zeros(n,n);
e(:,medoids) = uij;
end