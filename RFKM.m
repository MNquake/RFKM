function [F, obj_RFKM, iter, center, S] = RFKM(F, r, S, X, Noise, center, K)
    % X d*n 输入数据
    % F n*c 模糊矩阵
    % r 模糊矩阵系数
    % S n*c 噪音矩阵
    % Noise 噪音点个数
    % center d*c 中心点
    
    [~, n] = size(X);       % 取数据量
    c = size(F, 2);         % 取类别数          
    maxiter = 30;           % 最大迭代次数
    
    obj_RFKM = zeros(1, maxiter);   %存储目标函数的输出值
    center = findCenter(F, X, r, n, c, K, S);
 
    G = F.^r;                                           % 模糊矩阵 n*c
    dist = distfcm(X', center');                        % 距离矩阵 n*c
%     paint(X,c,F,S,center,1, dist);
    obj_RFKM(1) = sum(sum((dist.^2).*(G.*S)));          % 目标函数

    for iter = 2 : maxiter
        
        %% 更新模糊矩阵F
        for i = 1 : n
            tmp = distfcm(X(:,i)',center'); % 算出每个点与各个中心点的距离
            if any(tmp' == 0)               % 有零说明有重合
                col = find(tmp == 0);       % 找出重合的中心点
                F(i,:) = 0;                 % 说明一定属于该中心
                F(i,col) = 1;               % 直接将该点与该中心的隶属度设为1，其余设为零
                continue;
            end
            
            D = tmp.^(2/(1-r)) * ones(c,1);         % 更新模糊矩阵F公式中的分母
            F(i,:) = (tmp.^(2/(1-r)))./D;           % 更新模糊矩阵
        end
        
        %% 更新噪音矩阵S
        G = F.^r;
        dist = distfcm(X', center');
        S = ones(n,c);

        tmp = sort(dist, 2);              % 每个点距离各个中心点的距离进行升序排序
        [~, idx] = sortrows(-tmp,1:c);    % 每个点距离最近它最近的中心点的距离，对该距离进行降序排序,若有两个点距离最近它们最近的中心点的距离相同
        %-tmp为了降序排序                  % 则比较它们两个距离它们第二个距离最近的中心点，若仍相同，比较第三个，以此类推进行比较。
        S(idx(1:Noise),:) = 0;            % 取前Noise个作为噪音点

        %% 更新中心点矩阵C

        center = findCenter(F, X, r, n, c, K, S);


        %% 计算目标函数
        G = F.^r;
        dist = distfcm(X', center');                            % 距离矩阵 n*c
        obj_RFKM(iter) = sum(sum((dist.^2).*(G.*S)));           % 目标函数
%         paint(X,c,F,S,center,iter, dist);
        if (abs(obj_RFKM(iter)-obj_RFKM(iter-1)) < 10^-5)
            break;
        end

    end

end