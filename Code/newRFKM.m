function [F, obj_RFKM, iter, c_id] = newRFKM(F, r, X, Noise, K, S)

    [~, n] = size(X);       % 取数据量
    c = size(F, 2);         % 取类别数          
    maxiter = 30;           % 最大迭代次数
    d = distfcm(X',X');

    c_id = randperm(n,c);
    dist = d(:,c_id);       % 每个点与中心点的距离

    F = updateF(n,c,dist,r);
    
    G = F.^r;
    obj_RFKM(1) = sum(sum((dist.^2).*(G.*S)));

    for iter = 2 : maxiter
    
        %% 更新center
        G = F.^r;
        c_id = updateCenter(G,X,c_id,d,S,c,K);
        dist = d(:,c_id);

        %% 更新F
        F = updateF(n,c,dist,r);

        %% 更新噪音矩阵
        G = F.^r;
        dist = d(:,c_id);
        S = ones(n,c);

        tmp = sort(dist, 2);              % 每个点距离各个中心点的距离进行升序排序
        [~, idx] = sortrows(-tmp,1:c);    % 每个点距离最近它最近的中心点的距离，对该距离进行降序排序,若有两个点距离最近它们最近的中心点的距离相同
        %-tmp为了降序排序                  % 则比较它们两个距离它们第二个距离最近的中心点，若仍相同，比较第三个，以此类推进行比较。
        S(idx(1:Noise),:) = 0;            % 取前Noise个作为噪音点
    
        
        G = F.^r;
        dist = d(:,c_id);
        obj_RFKM(iter) = sum(sum((dist.^2).*(G.*S)));
%         paint(X,c,F,S,center,iter, dist);
        if (abs(obj_RFKM(iter)-obj_RFKM(iter-1)) < 10^-5)
            break;
        end
    end

end