function F = updateF(n,c,dist,r)
    tmp = dist.^(-1/(r-1));
    F = tmp./(sum(tmp,2)*ones(1,c));
    for i = 1 : n
        d = dist(i,:);
        if any(d' == 0)               % 有零说明有重合
            col = find(d == 0);       % 找出重合的中心点
            F(i,:) = 0;                 % 说明一定属于该中心
            F(i,col) = 1;               % 直接将该点与该中心的隶属度设为1，其余设为零
        end  
    end
end








% function F = updateF(n,c,dist,r)
%     F = zeros(n,c);
%     for i = 1:n
%         soma = 0;
%         for j = 1:c
%             soma = soma + dist(i,j)^(1/(r-1));
%         end
%         for j = 1:c
%             if dist(i,j) == 0
%                 F(i,j) = 1;
%             else
%                 F(i,j) = (dist(i,j)^(1/(r-1)))/soma;
%             end
%         end
%     end
% 
%     ff = zeros(n,c);
%     for i = 1 : n
%         d = dist(i,:);
%         if any(d' == 0)               % 有零说明有重合
%             col = find(d == 0);       % 找出重合的中心点
%             F(i,:) = 0;                 % 说明一定属于该中心
%             F(i,col) = 1;               % 直接将该点与该中心的隶属度设为1，其余设为零
%             continue;
%         end  
%         D = d.^(2/(1-r)) * ones(c,1);         % 更新模糊矩阵F公式中的分母
%         ff(i,:) = (d.^(2/(1-r)))./D;           % 更新模糊矩阵
%     end
%     t = 1;
%     ff = zeros(n,c);
%     for i = 1:n
%         soma = 0;
%         for j = 1:c
%             soma = soma + 1/dist(i,j)^(1/(r-1));
%         end
%         for j = 1:c
%             if dist(i,j) == 0
%                 ff(i,j) = 1;
%             else
%                 ff(i,j) = 1/(soma*dist(i,j)^(1/(r-1)));
%             end
%         end
%     end
% 
%     F_ = dist./(sum(dist,2)*ones(1,c));
% 
%     t = 1;
% end

