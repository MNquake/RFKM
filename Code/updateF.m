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

