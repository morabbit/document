tic
clear
%%均线策略
%%
%%导入数据
getfilename = ls('G:\量化投资策略编写\均线策略\数据\*.xlsx');
filename = cellstr(getfilename);
num = length(filename);
[data,str,all] = xlsread(filename{1});    %读取excel文件
stock{num}={filename{num},data,text};
Length = length(data);  %获取数据长度
Open = data(1:Length,1);    %获取开盘价数据
High = data(1:Length,2);    %获取最高价数据
Low = data(1:Length,3); %获取最低价数据
Close = data(1:Length,4);   %获取收盘价数据
%%
%%设置移动平均线
for x = 10;
    MA_PMT = x;    %设置移动平均线参数,PMT(parameter,参数的意思)
    MA = zeros(Length,1);  %给移动平均数据开辟空间
    [MA] = movavg(Close(:,1),MA_PMT,MA_PMT);%计算移动平均值
    MA(1:MA_PMT-1,1) = Close(1:MA_PMT-1,1); %将错误值更改为相应时间的收盘价
    Date = str(3:end,1);    %读取日期
    %%
    %%持仓点位的设置及市值计算
    Pos = zeros(Length,1);  %给仓位开辟空间
    %设置开多及平仓条件及仓位填充,问题：储存开多的下标，计算胜率用
    for t = MA_PMT-1:Length
        if Close(t-1) > MA(t-1)
            Pos(t) = 1; %满足开多条件
            continue
        end
        if Close(t-1) <= MA(t-1)
           Pos(t) = 0;  %满足平仓条件
           continue
        end
        if Pos(t-1) == 1
            Pos(t) = 1; %多头仓位填充
        elseif Pos(t-1) == 0
            Pos(t) = 0; %平仓仓位填充
        end
    end
    Pos(Length) = 0; %最后一个时点平仓

    %买卖点位的设置
    Time = zeros(Length,1); %先给买卖点位开辟空间
    for i = 2:Length
        Time(i) = Pos(i) - Pos(i-1);    %判断开多和平仓时点
    end

    %生成以最高价买以最低价卖的市值序列及手续费
    MV = zeros(Length,1);   %给市值序列开辟空间
    MV(1,1) = 100000;   %初始投入10万元本金
    Trans_Cost = 0;   %给交易手续费开辟空间
    yj = 0.0002;    %佣金
    ghf = 0.00002;    %过户费
    yhs = 0.001;    %印花税
    for i = 2:Length
        if Pos(i) == 0 && Time(i) == 0
            MV(i) = MV(i-1);
        elseif Pos(i) == 1 && Time(i) == 1
            MV(i) = MV(i-1)/High(i)*Close(i);
            Trans_Cost = MV(i)*(yj+ghf);
        elseif Pos(i) == 1 && Time(i) == 0
            MV(i) = MV(i-1)/Close(i-1)*Close(i);
        else
            MV(i) = MV(i-1)/Close(i-1)*Low(i)*(1-yj-ghf-yhs)-Trans_Cost;
        end
    end
    %%
    %%总收益率、年化收益率、每年收益率、收益标准差、
    ...最大回撤、最长回撤期，盈亏比，胜率，夏普比率；

    %计算总收益率
    TalRet = MV(end,1)/MV(1)-1;

    %年化收益率
    Period = datenum(Date(end))-datenum(Date(1));
    AvgRet = (1+TalRet)^(365/Period)-1;

    %每年收益率
    Year = zeros(Length,1);
    Year(1,1) = 1;
    Year(end,1) = 1;
    for i = 3:Length-1
        Year(i-1) = year(Date(i)) - year(Date(i-1));
    end
    Loc = find(Year==1);
    for i = Loc
        RetPY = MV(i);
    end
    PY = zeros(length(RetPY),1);
    PY = num2cell(PY);
    for i = 1:length(RetPY)
        PY(i) =  Date(Loc(i));
    end
    RoRPY = zeros(length(RetPY),1);
    for i = 2:length(RetPY)
        RoRPY(i) = log(RetPY(i)/RetPY(i-1));
    end

    %每年收益标准差
    Std_RoRPY = std(RoRPY);

    %最大回撤
    MaxDraw = zeros(Length,1);
    for i = 2:Length
        MaxDraw(i) = (max(MV(1:i,1))-MV(i))/max(MV(1:i,1));
    end
    MaxDraw_Value = max(MaxDraw);

    %最长回撤期
    MaxDrawDuration = zeros(Length,1);
    for i = 2:Length
        if MV(i)<max(MV(1:i,1)) == 1
            MaxDrawDuration(i) = MaxDrawDuration(i-1)+1;
        else
            MaxDrawDuration(i) = 0;
        end
    end
    MaxDrawDuration_Value = max(MaxDrawDuration);

    %计算胜率
    Location = find(Pos == 1 & Time ==1 | Pos == 0 & Time == -1);
    Begin = zeros(length(Location)/2,1);
    End = zeros(length(Location)/2,1);
    for k = 1:length(Location)
        if mod(k,2) ~= 0
            j = (k+1)/2;
            Begin(j) = MV(Location(k));
        else
            j = k/2;
            End(j) = MV(Location(k));
        end
    end
    Profit = End-Begin;
    Positive_Return = length(find(Profit>0));
    Negative_Return = length(find(Profit<0));
    Win_Percentage = Positive_Return/(Negative_Return + ...
        Positive_Return);

    %计算盈亏比
    Positive_ReturnSum = sum(Profit(Profit>0));
    Negative_ReturnSum = sum(Profit(Profit<0));
    Profit_to_Loss = Positive_ReturnSum/abs(Negative_ReturnSum);

    %计算夏普比率
    RoR = (MV(2:end) - MV(1:end-1))./MV(1:end-1);
    SR = mean(RoR)/std(RoR);

    %%
    %将结果输出值excel
    filename = 'Results.xlsx';
    A = {'总收益率',TalRet;'年化收益率',AvgRet;...
        '最大回撤',MaxDraw_Value;'最长回撤期',MaxDrawDuration_Value;...
        '胜率',Win_Percentage;'盈亏比',Profit_to_Loss;'夏普比率',SR};
    sheet = x - 9;
    xlrange = 'A2';
    xlswrite(filename,A,sheet,xlrange);
    xlswrite(filename,PY,sheet,'A9');
    xlswrite(filename,RoRPY,sheet,'B9');
    xlswrite(filename,MA_PMT,sheet,'B1');
end
toc