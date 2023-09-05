function c_id = updateCenter(G,X,c_id,d,S,c,K)
    
    center = ((G.*S)'*X')./(sum((G.*S)',2)*ones(1,size(X',2)));
    center = center';

    dist = distfcm(X', center');         % k-means的中心距离矩阵
    
    [~, col_idx] = max(G, [], 2);

    for j = 1 : c
    
        [row_idx, ~] = find(col_idx == j);          % row_idx就记录了哪些点属于第j类
        Attr_dist = dist(row_idx',j);               % 取出属于这一类的点与中心点的距离
        if size(Attr_dist,1) == 0                   % 初始化的时候，F矩阵有可能会让某一类没有点   
            continue;
        end
        [~, idx] = sort(Attr_dist);                 % 将距离进行排序，idx记录了原本在row_idx中的位置
        [col, ~] = size(Attr_dist);
        if col > K                                  % 要判断一下属于该中心点的样本量是否大于遍历的个数
            len = K;
        else
            len = col;
        end
                         % row_idx(idx(1)): idx存放的是该点原本在row_idx中的位置，而row_idx记录的是第几个点
        t = 1;
        while S(row_idx(idx(t)),1) == 0                     % 处理噪音点 直接跳过噪音点       
            t = t + 1;
            if t > size(idx,1)
                break;
            end
        end
        if t > size(idx,1)
            continue;
        end
        
        c_id(j) = row_idx(idx(t));
        c_id_ = c_id;
        dd = d(:,c_id_);       
        minres = sum(sum((dd.^2).*G)); % 第一个点的目标函数输出值
%         minindices = t;
        for i = t + 1 : len
            if S(row_idx(idx(i)),1) == 0                    % 处理噪音点 直接跳过噪音点       
                continue;
            end
            c_id_(j) = row_idx(idx(i));
            dd = d(:,c_id_);
            if minres > sum(sum((dd.^2).*G))       
                minres = sum(sum((dd.^2).*G)); % 有更小的就更新
%                 minindices = i;
                c_id(j) = c_id_(j);
            end
        end
    end

end
