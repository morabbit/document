tic
clear
%%���߲���
%%
%%��������
[data,str,all] = xlsread('G:\����Ͷ�ʲ��Ա�д\���߲���\����\159919.xlsx');
Length = length(data);  %��ȡ���ݳ���
Open = data(1:Length,1);    %��ȡ���̼�����
High = data(1:Length,2);    %��ȡ��߼�����
Low = data(1:Length,3); %��ȡ��ͼ�����
Close = data(1:Length,4);   %��ȡ���̼�����
Names = {'��������','�껯������','ÿ�������׼��',...   %Ϊ���ݵ�������׼��
    '���س�','��س���','ʤ��',...
    'ӯ����','���ձ���'};
T = table('RowName',Names);
%%
%%�����ƶ�ƽ���ߣ����е�һ�β���Ѱ��
for x = [3,10:10:300]
    MA_PMT = x;    %�����ƶ�ƽ���߲���,PMT(parameter,��������˼)
    MA = zeros(Length,1);  %���ƶ�ƽ�����ݿ��ٿռ�
    [MA] = movavg(Close(:,1),MA_PMT,MA_PMT);%�����ƶ�ƽ��ֵ
    MA(1:MA_PMT-1,1) = Close(1:MA_PMT-1,1); %������ֵ����Ϊ��Ӧʱ������̼�
    Date = str(3:end,1);    %��ȡ����
    %��ͼ
    % plot(Close,'r');
    % hold on;
    % plot(MA,'g:')
    % set(gca,'xtick',1:300:Length);
    % set(gca,'xticklabel',Date(1:300:Length));
    % legend('Close','MA');
    % title('�ƶ�ƽ������');
    %%
    %%�ֲֵ�λ�����ü���ֵ����
    Pos = zeros(Length,1);  %����λ���ٿռ�
    %���ÿ��༰ƽ����������λ���,���⣺���濪����±꣬����ʤ����
    for t = MA_PMT-1:Length
        if Close(t-1) > MA(t-1)
            Pos(t) = 1; %���㿪������
            continue
        end
        if Close(t-1) <= MA(t-1)
           Pos(t) = 0;  %����ƽ������
           continue
        end
        if Pos(t-1) == 1
            Pos(t) = 1; %��ͷ��λ���
        elseif Pos(t-1) == 0
            Pos(t) = 0; %ƽ�ֲ�λ���
        end
    end
    Pos(Length) = 0; %���һ��ʱ��ƽ��

    %������λ������
    Time = zeros(Length,1); %�ȸ�������λ���ٿռ�
    for i = 2:Length
        Time(i) = Pos(i) - Pos(i-1);    %�жϿ����ƽ��ʱ��
        if Time(i) == 1
    %         plot(i,Close(i),'b*','markersize',5);   %��ͼ�б�ע����ʱ��
        elseif Time(i) == -1
    %         plot(i,Close(i),'ko','markersize',5);   %��ͼ�б�עƽ��ʱ��
        end
    end

    %��������߼�������ͼ�������ֵ���м�������
    MV = zeros(Length,1);   %����ֵ���п��ٿռ�
    MV(1,1) = 100000;   %��ʼͶ��10��Ԫ����
    Trans_Cost = 0;   %�����������ѿ��ٿռ�
    yj = 0.0002;    %Ӷ��
    ghf = 0.00002;    %������
    yhs = 0.001;    %ӡ��˰
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
    %%�������ʡ��껯�����ʡ�ÿ�������ʡ������׼�
    ...���س�����س��ڣ�ӯ���ȣ�ʤ�ʣ����ձ��ʣ�

    %������������
    TalRet = MV(end,1)/MV(1)-1;

    %�껯������
    Period = datenum(Date(end))-datenum(Date(1));
    AvgRet = (1+TalRet)^(365/Period)-1;

    %ÿ��������
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

    %ÿ�������׼��
    Std_RoRPY = std(RoRPY);

    %���س�
    MaxDraw = zeros(Length,1);
    for i = 2:Length
        MaxDraw(i) = (max(MV(1:i,1))-MV(i))/max(MV(1:i,1));
    end
    MaxDraw_Value = max(MaxDraw);

    %��س���
    MaxDrawDuration = zeros(Length,1);
    for i = 2:Length
        if MV(i)<max(MV(1:i,1)) == 1
            MaxDrawDuration(i) = MaxDrawDuration(i-1)+1;
        else
            MaxDrawDuration(i) = 0;
        end
    end
    MaxDrawDuration_Value = max(MaxDrawDuration);

    %����ʤ��
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

    %����ӯ����
    Positive_ReturnSum = sum(Profit(Profit>0));
    Negative_ReturnSum = sum(Profit(Profit<0));
    Profit_to_Loss = Positive_ReturnSum/abs(Negative_ReturnSum);

    %�������ձ���
    RoR = (MV(2:end) - MV(1:end-1))./MV(1:end-1);
    SR = mean(RoR)/std(RoR);

    %%
    %��������table
    eval(['MA' num2str(MA_PMT) '=' 'zeros' '(' num2str(8) ',' num2str(1) ')' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(1) ',' num2str(1) ')' '=' 'TalRet' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(2) ',' num2str(1) ')' '=' 'AvgRet' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(3) ',' num2str(1) ')' '=' 'Std_RoRPY' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(4) ',' num2str(1) ')' '=' 'MaxDraw_Value' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(5) ',' num2str(1) ')' '=' 'MaxDrawDuration_Value' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(6) ',' num2str(1) ')' '=' 'Win_Percentage' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(7) ',' num2str(1) ')' '=' 'Profit_to_Loss' ';'])
    eval(['MA' num2str(MA_PMT) '(' num2str(8) ',' num2str(1) ')' '=' 'SR' ';'])
    eval(['MA' num2str(MA_PMT) '=' 'table' '(' 'MA' num2str(MA_PMT) ')']);
    T = [T,eval(['MA' num2str(MA_PMT)])];
end
[C, I] = max(T{1,:});
optimal_MA_PMT = cell2mat(T.Properties.VariableNames(I))
optimal_MA_PMT_Number = str2num(cell2mat(regexp(optimal_MA_PMT, '\d', 'match')))

%%
%%���ж��β���Ѱ��
T = table('RowName',Names);
if optimal_MA_PMT_Number-10 < 0
    for x = [3:1:optimal_MA_PMT_Number+10]
        MA_PMT = x;    %�����ƶ�ƽ���߲���,PMT(parameter,��������˼)
        MA = zeros(Length,1);  %���ƶ�ƽ�����ݿ��ٿռ�
        [MA] = movavg(Close(:,1),MA_PMT,MA_PMT);%�����ƶ�ƽ��ֵ
        MA(1:MA_PMT-1,1) = Close(1:MA_PMT-1,1); %������ֵ����Ϊ��Ӧʱ������̼�
        Date = str(3:end,1);    %��ȡ����
        %��ͼ
        % plot(Close,'r');
        % hold on;
        % plot(MA,'g:')
        % set(gca,'xtick',1:300:Length);
        % set(gca,'xticklabel',Date(1:300:Length));
        % legend('Close','MA');
        % title('�ƶ�ƽ������');
        %%
        %%�ֲֵ�λ�����ü���ֵ����
        Pos = zeros(Length,1);  %����λ���ٿռ�
        %���ÿ��༰ƽ����������λ���,���⣺���濪����±꣬����ʤ����
        for t = MA_PMT-1:Length
            if Close(t-1) > MA(t-1)
                Pos(t) = 1; %���㿪������
                continue
            end
            if Close(t-1) <= MA(t-1)
               Pos(t) = 0;  %����ƽ������
               continue
            end
            if Pos(t-1) == 1
                Pos(t) = 1; %��ͷ��λ���
            elseif Pos(t-1) == 0
                Pos(t) = 0; %ƽ�ֲ�λ���
            end
        end
        Pos(Length) = 0; %���һ��ʱ��ƽ��

        %������λ������
        Time = zeros(Length,1); %�ȸ�������λ���ٿռ�
        for i = 2:Length
            Time(i) = Pos(i) - Pos(i-1);    %�жϿ����ƽ��ʱ��
            if Time(i) == 1
        %         plot(i,Close(i),'b*','markersize',5);   %��ͼ�б�ע����ʱ��
            elseif Time(i) == -1
        %         plot(i,Close(i),'ko','markersize',5);   %��ͼ�б�עƽ��ʱ��
            end
        end

        %��������߼�������ͼ�������ֵ���м�������
        MV = zeros(Length,1);   %����ֵ���п��ٿռ�
        MV(1,1) = 100000;   %��ʼͶ��10��Ԫ����
        Trans_Cost = 0;   %�����������ѿ��ٿռ�
        yj = 0.0002;    %Ӷ��
        ghf = 0.00002;    %������
        yhs = 0.001;    %ӡ��˰
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
        %%�������ʡ��껯�����ʡ�ÿ�������ʡ������׼�
        ...���س�����س��ڣ�ӯ���ȣ�ʤ�ʣ����ձ��ʣ�

        %������������
        TalRet = MV(end,1)/MV(1)-1;

        %�껯������
        Period = datenum(Date(end))-datenum(Date(1));
        AvgRet = (1+TalRet)^(365/Period)-1;

        %ÿ��������
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

        %ÿ�������׼��
        Std_RoRPY = std(RoRPY);

        %���س�
        MaxDraw = zeros(Length,1);
        for i = 2:Length
            MaxDraw(i) = (max(MV(1:i,1))-MV(i))/max(MV(1:i,1));
        end
        MaxDraw_Value = max(MaxDraw);

        %��س���
        MaxDrawDuration = zeros(Length,1);
        for i = 2:Length
            if MV(i)<max(MV(1:i,1)) == 1
                MaxDrawDuration(i) = MaxDrawDuration(i-1)+1;
            else
                MaxDrawDuration(i) = 0;
            end
        end
        MaxDrawDuration_Value = max(MaxDrawDuration);

        %����ʤ��
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

        %����ӯ����
        Positive_ReturnSum = sum(Profit(Profit>0));
        Negative_ReturnSum = sum(Profit(Profit<0));
        Profit_to_Loss = Positive_ReturnSum/abs(Negative_ReturnSum);

        %�������ձ���
        RoR = (MV(2:end) - MV(1:end-1))./MV(1:end-1);
        SR = mean(RoR)/std(RoR);

        %%
        %��������table
        eval(['MA' num2str(MA_PMT) '=' 'zeros' '(' num2str(8) ',' num2str(1) ')' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(1) ',' num2str(1) ')' '=' 'TalRet' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(2) ',' num2str(1) ')' '=' 'AvgRet' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(3) ',' num2str(1) ')' '=' 'Std_RoRPY' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(4) ',' num2str(1) ')' '=' 'MaxDraw_Value' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(5) ',' num2str(1) ')' '=' 'MaxDrawDuration_Value' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(6) ',' num2str(1) ')' '=' 'Win_Percentage' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(7) ',' num2str(1) ')' '=' 'Profit_to_Loss' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(8) ',' num2str(1) ')' '=' 'SR' ';'])
        eval(['MA' num2str(MA_PMT) '=' 'table' '(' 'MA' num2str(MA_PMT) ')']);
        T = [T,eval(['MA' num2str(MA_PMT)])];
    end
    [C, I] = max(T{1,:});
    optimal_MA_PMT = cell2mat(T.Properties.VariableNames(I))
    optimal_MA_PMT_Number = str2num(cell2mat(regexp(optimal_MA_PMT, '\d', 'match')))    %���Ų���Ҫ���
    Daily_Return = zeros(Length,1);
    Daily_Return = (Close(2:end)-Close(1:end-1))./Close(1:end-1)
    Var_Daily_Return = var(Daily_Return) %����Ҫ���
else
    for x = [optimal_MA_PMT_Number-10:1:optimal_MA_PMT_Number+10]
       MA_PMT = x;    %�����ƶ�ƽ���߲���,PMT(parameter,��������˼)
        MA = zeros(Length,1);  %���ƶ�ƽ�����ݿ��ٿռ�
        [MA] = movavg(Close(:,1),MA_PMT,MA_PMT);%�����ƶ�ƽ��ֵ
        MA(1:MA_PMT-1,1) = Close(1:MA_PMT-1,1); %������ֵ����Ϊ��Ӧʱ������̼�
        Date = str(3:end,1);    %��ȡ����
        %��ͼ
        % plot(Close,'r');
        % hold on;
        % plot(MA,'g:')
        % set(gca,'xtick',1:300:Length);
        % set(gca,'xticklabel',Date(1:300:Length));
        % legend('Close','MA');
        % title('�ƶ�ƽ������');
        %%
        %%�ֲֵ�λ�����ü���ֵ����
        Pos = zeros(Length,1);  %����λ���ٿռ�
        %���ÿ��༰ƽ����������λ���,���⣺���濪����±꣬����ʤ����
        for t = MA_PMT-1:Length
            if Close(t-1) > MA(t-1)
                Pos(t) = 1; %���㿪������
                continue
            end
            if Close(t-1) <= MA(t-1)
               Pos(t) = 0;  %����ƽ������
               continue
            end
            if Pos(t-1) == 1
                Pos(t) = 1; %��ͷ��λ���
            elseif Pos(t-1) == 0
                Pos(t) = 0; %ƽ�ֲ�λ���
            end
        end
        Pos(Length) = 0; %���һ��ʱ��ƽ��

        %������λ������
        Time = zeros(Length,1); %�ȸ�������λ���ٿռ�
        for i = 2:Length
            Time(i) = Pos(i) - Pos(i-1);    %�жϿ����ƽ��ʱ��
            if Time(i) == 1
        %         plot(i,Close(i),'b*','markersize',5);   %��ͼ�б�ע����ʱ��
            elseif Time(i) == -1
        %         plot(i,Close(i),'ko','markersize',5);   %��ͼ�б�עƽ��ʱ��
            end
        end

        %��������߼�������ͼ�������ֵ���м�������
        MV = zeros(Length,1);   %����ֵ���п��ٿռ�
        MV(1,1) = 100000;   %��ʼͶ��10��Ԫ����
        Trans_Cost = 0;   %�����������ѿ��ٿռ�
        yj = 0.0002;    %Ӷ��
        ghf = 0.00002;    %������
        yhs = 0.001;    %ӡ��˰
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
        %%�������ʡ��껯�����ʡ�ÿ�������ʡ������׼�
        ...���س�����س��ڣ�ӯ���ȣ�ʤ�ʣ����ձ��ʣ�

        %������������
        TalRet = MV(end,1)/MV(1)-1;

        %�껯������
        Period = datenum(Date(end))-datenum(Date(1));
        AvgRet = (1+TalRet)^(365/Period)-1;

        %ÿ��������
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

        %ÿ�������׼��
        Std_RoRPY = std(RoRPY);

        %���س�
        MaxDraw = zeros(Length,1);
        for i = 2:Length
            MaxDraw(i) = (max(MV(1:i,1))-MV(i))/max(MV(1:i,1));
        end
        MaxDraw_Value = max(MaxDraw);

        %��س���
        MaxDrawDuration = zeros(Length,1);
        for i = 2:Length
            if MV(i)<max(MV(1:i,1)) == 1
                MaxDrawDuration(i) = MaxDrawDuration(i-1)+1;
            else
                MaxDrawDuration(i) = 0;
            end
        end
        MaxDrawDuration_Value = max(MaxDrawDuration);

        %����ʤ��
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

        %����ӯ����
        Positive_ReturnSum = sum(Profit(Profit>0));
        Negative_ReturnSum = sum(Profit(Profit<0));
        Profit_to_Loss = Positive_ReturnSum/abs(Negative_ReturnSum);

        %�������ձ���
        RoR = (MV(2:end) - MV(1:end-1))./MV(1:end-1);
        SR = mean(RoR)/std(RoR);

        %%
        %��������table
        eval(['MA' num2str(MA_PMT) '=' 'zeros' '(' num2str(8) ',' num2str(1) ')' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(1) ',' num2str(1) ')' '=' 'TalRet' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(2) ',' num2str(1) ')' '=' 'AvgRet' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(3) ',' num2str(1) ')' '=' 'Std_RoRPY' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(4) ',' num2str(1) ')' '=' 'MaxDraw_Value' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(5) ',' num2str(1) ')' '=' 'MaxDrawDuration_Value' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(6) ',' num2str(1) ')' '=' 'Win_Percentage' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(7) ',' num2str(1) ')' '=' 'Profit_to_Loss' ';'])
        eval(['MA' num2str(MA_PMT) '(' num2str(8) ',' num2str(1) ')' '=' 'SR' ';'])
        eval(['MA' num2str(MA_PMT) '=' 'table' '(' 'MA' num2str(MA_PMT) ')']);
        T = [T,eval(['MA' num2str(MA_PMT)])];
    end 
    [C, I] = max(T{1,:});
    optimal_MA_PMT = cell2mat(T.Properties.VariableNames(I))
    optimal_MA_PMT_Number = str2num(cell2mat(regexp(optimal_MA_PMT, '\d', 'match')))    %���Ų���Ҫ���
    Daily_Return = zeros(Length,1);
    Daily_Return = (Close(2:end,1)-Close(1:end-1,1))./Close(1:end-1,1)
    Var_Daily_Return = var(Daily_Return) %����Ҫ���
end
%%
%%��������ŵľ��߲������������ʷ������
Final_Results_RowName = {'���о��߲���'; '����'};
Final_Results = table('RowName', Final_Results_RowName)
b = cell2mat(str(1,1));
Asset_Name = str2num(cell2mat(regexp(b, '\d', 'match')));
eval(['Asset' num2str(Asset_Name) '(' num2str(1) ',' num2str(1) ')' '=' 'optimal_MA_PMT_Number'])
eval(['Asset' num2str(Asset_Name) '(' num2str(2) ',' num2str(1) ')' '=' 'Var_Daily_Return'])
eval(['Result' '=' 'table' '(' 'Asset' num2str(Asset_Name) ')']);
Final_Results = [Final_Results, Result]
%��������Ϊ�������
toc