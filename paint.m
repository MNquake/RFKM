function paint(data, c, F, S, center, iter, dist)

    [n, ~] = size(S);
    
    for i = 1:n
        [~, idx] = sort(F(i,:), 'descend');         % 将每一行从大到小排序
        tmp(i,:) = idx;                             % 记录原始的列号，也就是这列数据的概率是哪一个中心点的
    end
    col_idx = tmp(:,1);                             % 第一列即记录了每个点属于哪一个中心点
    hold on;
    color = ['r','b','g','k'];
    for j = 1 : c
        [row_idx, ~] = find(col_idx == j);          % row_idx就记录了哪些点属于第j类
        scatter(data(1,row_idx),data(2,row_idx),'o',color(j));
        scatter(center(1,j),center(2,j),'fill',color(j));
    end

    S_ = S(:,1);
    rows = find(S_ == 0);
    scatter(data(1,rows),data(2,rows), 'v','y', 'LineWidth',1.2);

    str = '迭代次数: ';
    title([str, num2str(iter)]);                    %标题
    axis([-23, 23, -23, 23]);
    grid on;
    
    filepath = 'G:\Code\RFKM\2\';
    
    filename = [num2str(iter), '.png'];
    saveas(gcf,[filepath, filename]);

    filename = ['dist-', num2str(iter), '.mat'];
    save([filepath, filename],'dist');

    filename = ['F-', num2str(iter), '.mat'];
    save([filepath, filename],'dist');
    
    clf;
end