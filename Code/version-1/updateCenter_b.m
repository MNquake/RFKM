function [c_id,F] = updateCenter(X,K,F,c_id,d,S,r,n,c)

    G = F.^r;
    center = (G'*X')./(sum(G',2)*ones(1,size(X',2)));
    center = center';                                         % k-means的聚类中心
    
    [~, col_idx] = max(F, [], 2);                             % col_idx即记录了每个点属于哪一个中心点 n*1
    dist = distfcm(X', center');                              % 各点与kmeans聚类中心的距离
    for j = 1:c
        [row_idx, ~] = find(col_idx == j);                  % row_idx记录了哪些点属于第j类
        Attr_dist = dist(row_idx',j);                       % 取出属于这一类的点与中心点的距离
        if size(Attr_dist,1) == 0                           % 初始化的时候，F矩阵有可能会让某一类没有点   
            continue;
        end
        [~, idx] = sort(Attr_dist);                 % 对距离进行排序，idx记录了原本在row_idx中的位置
        [col, ~] = size(Attr_dist);
        if col > K                                  % 要判断一下属于该中心点的样本量是否大于遍历的个数
            len = K;
        else
            len = col;
        end
                   
        % row_idx(idx(1)): idx存放的是该点原本在row_idx中的位置，而row_idx记录的是第几个点
        md = dist(c_id(j), j);    % 先把自己点作为最小点，认为是中心点，用于遍历时的比较        
        cid_ = c_id;              % 暂时记录
        cid_(j) = c_id(j);
        F_ = updateF(n,c,d(:,cid_),r);               %在查找点的过程中更新F
        G = F.^r;
        minres = sum(sum((md.^2).*(G(c_id(j),j)))); % 第一个点的目标函数输出值
        minindices = -1;
        for i = 1 : len
            if S(row_idx(idx(i)),1) == 0                    % 处理噪音点 直接跳过噪音点       
                continue;
            end
            md = dist(row_idx(idx(i)), j);
            cid_ = c_id;
            cid_(j) = row_idx(idx(i));
            F_ = updateF(n,c,d(:,cid_),r);
            G = F_.^r;
            rr = sum(sum((md.^2).*(G(row_idx(idx(i)),j))));
            if minres > rr        
                minres = rr; % 有更小的就更新
                minindices = i;
                F = F_;
            end
        end
        if minindices ~= -1
            c_id(j) = row_idx(idx(minindices));
        end
    end
    
end