%测试并行计算的原理
clear;
clc;
myCluster = parcluster('local');
myCluster.NumWorkers = 4;  % 'Modified' property now TRUE
saveProfile(myCluster);    % 'local' profile now updated,
                           % 'Modified' property now FALSE
para=parpool('local',4);



miu=0.9;
delta=0.1; 
beta=0.99;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
iter=1;
l_k=0;%以上是设置参数
[f_name,Pnameh]=load_file();
load(f_name);
p_num=1;
accuracy=0.000001;
%link_capacity=link_capacity;
t_alfa=0.15;%BPR函数的系数
t_beta=4;%BPR函数的指数
w=length(OD(:,1));%OD对的数量
link_num=length(link(:,1));%link的数量
al=0;
be=0;
ac=0;
acc=0;
p_size=p_num*ones(w,1,'double');%每个OD对有多少条路径
%load('CGsf2.mat');%加载路网，内含各种输入矩阵-OD、link、route_link、OD_route
par_t=[ 1;0.8;0.8;0.4;0.4];%time value参数
par_i=[ 0;0.3;0.4;0.3;0.4];%inconvenience参数
par_p=[ 0;-1;-1;1;1];%price参数
bench=20;%基准价格
par_c=[ 0;0.5;0.5;0;0];%RD补偿参数
par_r=[ 0;0;0;0.1;0.1];%R价格参数
par_fix=[1;1;1;0;0];%trip费用
r1 = [0,1,0,-1,0];%R1的匹配约束
r2 = [0,0,1,0,-0.5];%R2的匹配约束
veh = [1,1,1,0,0];%车辆的求和矩阵


%link_number=size(link,1);%得到路段数
%OD_number=size(OD,1);%得到OD数
%OD_demand=zeros(OD_number,1,'double');%OD需求量
%OD_route_number=zeros(OD_number,1,'double');%每个OD对应的路径数
%link_volume=zeros(link_number,1,'double');%路段流量
%link_volume_pre=zeros(link_number,1,'double');%前一次的路段流量
%link_cost=zeros(link_number,1,'double');%路段费用
%d_link_cost=zeros(link_number,1,'double');%路段费用的导数
%link_capacity=zeros(link_number,1,'double');%路段通行能力
%link_t0=zeros(link_number,1,'double');%路段自由流旅行时间
%link_b0=zeros(link_number,1,'double');%路段容量约束对应的拉格朗日乘子
%link_r0=zeros(link_number,1,'double');%增强拉格朗日乘子法中的惩罚因子
%link_origin=zeros(link_number,1,'double');%路段的起点
%link_destination=zeros(link_number,1,'double');%路段的终点
%link_zero=zeros(link_number,1,'double');%用于计算拉格朗日乘子
%link_Ca_cost=zeros(link_number,1,'double');%路段在达到通行能力时的费用
%route_link_1=zeros(1,link_number,'double');%路段在达到通行能力时的费用
%x_change=zeros(link_number,1,'double');%路段流量每次的变化量
%min_route=zeros(OD_number,1,'double');%存储每次迭代计算出的各个OD对的最短路径
%min_route_cost=zeros(OD_number,1,'double');%存储每次迭代计算出的各个OD对的最短路径

%以下更新初始路径费用


    for i=1:w %
        B_c{i}=reshape(pflow{i},5,[]);%把每条测试的路径流量按列排列
        veh_pflow_c{i}=veh*B_c{i};%计算测试的车辆路径流量
        x_c(:,i)=path{i}'*veh_pflow_c{i}';%计算测试的层流             
        ODflow{i}=sum(B_c{i},2);%计算测试的OD流
    end    
    lflow=sum(x_c,2);  %计算测试的路段流  
    lcost=link_t0.*(1+t_alfa*(lflow./link_capacity).^t_beta);%计算测试的路段费用
    for i=1:w %
        tpc{i}=(par_t+par_i)*(path{i}*lcost)';%计算测试的路径费用的拥堵费用部分
        mpc{i}=repmat(par_p.*(bench+(par_r-par_c).*ODflow{i}),p_size(i),1);%计算测试的路径费用的货币费用部分
        pcost{i}=tpc{i}(:)+mpc{i}+repmat(par_fix,p_size(i),1);%计算测试的路径费用
    end

%load('column1_0328_2.5OD_1304iter_Rgap0.001_completed.mat');%断点续传：把上面的都%掉
accuracy=0.01;

    
for iter=1:100000
    tic;
%求最短路
DG=sparse(link_origin',link_destination',lcost);
for i=1:length(DG)
    for j=1:length(DG)
        if DG(i,j)==0
            DG(i,j)=Inf;
        end
    end
end
route_num=1;%路径记号
route_link=zeros(w,link_num);
for i=1:w
[routeCell,Dist]=kShortestPath(DG,OD(i,1),OD(i,2),1);%最短路算法，求出来的是[1条最短路径（以结点为标记，cell）,费用]
    i_path_num=1;
    for ci =1:length(routeCell)
        route=cell2mat(routeCell(ci));%把第ci条路径的结点表示成mat
        for j=1:size(route,2)-1
            [link_ID,l]=find(link(:,5)==route(j)*100+route(j+1));%把以路段为标记的最短路转换成以关联矩阵为标记
            route_link(route_num,link_ID)=1;%生成关联矩阵
        end        
        OD_route(i,route_num)=1;
        route_num=route_num+1;        
    end
    %OD_route_number(i)=route_num;
end
%最短路求完了
%route_link(1,2)=1;%测试，把OD对1的最短路改掉

%下面比对路径集中是否有最短路，如果有就忽略，如果没有就把最短路添加进路径集
for i=1:w
    not_match=0;%重置路径不同的计数
    for p=1:p_size(i) %在当前OD对的路径集里面进行比对
       lib=find(path{i}(p,:)==route_link(i,:)); %判断当前路径是否与最短路相同
       if size(lib,2)==link_num
           break;%如果从路径集里找到了最短路就停止比对，继续比对下一个OD对
       else
           not_match=not_match+1;%如果当前路径不同的话，路径不同的计数加1
       end 
    end 
    if not_match==p_size(i) %如果最短路不在当前路径集里
        p_size(i)=p_size(i)+1;%就路径数量加1
        path{i}(p_size(i),:)=route_link(i,:);%把最短路放进路径集  
        pflow{i}(p_size(i)*5-4:p_size(i)*5,1)=0;%把最短路的路径流量（0流量，5种用户）放进路径流量集
        pcost{i}(p_size(i)*5-4:p_size(i)*5,1)=(par_t+par_i).*(path{i}(p_size(i),:)*lcost)+par_p.*(bench+(par_r-par_c).*ODflow{i})+par_fix;%把最短路的路径费用放进费用集里
        %下面更新约束
        H{i}=eye(5*p_size(i));
        p_beq{i}=zeros(2*p_size(i),1);
        beq{i}=[OD_demand(i);p_beq{i}];
        R1(i,1:p_size(i)) = {r1};
        R2(i,1:p_size(i)) = {r2};
        A1{i} = blkdiag(R1{i,:});%生成流量约束1的系数矩阵（p_size维分块对角矩阵，每个分块等于[0,1,0,-1,0]）
        A2{i} = blkdiag(R2{i,:});%生成流量约束2的系数矩阵（p_size维分块对角矩阵，每个分块等于[0,0,1,0,-0.5]）
        d_Aeq{i}=ones(1,5*p_size(i));%需求约束的系数
        Aeq{i}=[d_Aeq{i};A1{i};A2{i};];%等式约束的系数
        lb{i}=zeros(5*p_size(i),1);%路径流量的下界
    end    
end


%判断精度
parfor i=1:w %
    pflow_e{i}=quadprog(H{i},1*pcost{i}-pflow{i},[],[],Aeq{i},beq{i},lb{i},[],pflow{i});%用二次规划来计算投影
    e(i)=norm(pflow{i}-pflow_e{i},2);
    f(i)=norm(pflow{i});
end

%下面计算精度并判断

ac=(e*e')^0.5/(f*f')^0.5;
acc(iter)=ac;%把第精度存到acc里
if ac<accuracy
    break;
end
%判断完毕


%下面找步长
for l_l=0:100
    beta_l=miu^l_l*beta;
    parfor i=1:w %
        pflow_l{i}=quadprog(H{i},beta_l*pcost{i}-pflow{i},[],[],Aeq{i},beq{i},lb{i},[],pflow{i});%用二次规划来计算投影
        B_l{i}=reshape(pflow_l{i},5,[]);%把每条测试的路径流量按列排列
        veh_pflow_l{i}=veh*B_l{i};%计算测试的车辆路径流量
        x_l(:,i)=path{i}'*veh_pflow_l{i}';%计算测试的层流             
        ODflow_l{i}=sum(B_l{i},2);%计算测试的OD流
    end
    lflow_l=sum(x_l,2);  %计算测试的路段流  
    lcost_l=link_t0.*(1+t_alfa*(lflow_l./link_capacity).^t_beta);%计算测试的路段费用
    
% 下面计算各OD对的步长判据    
    for i=1:w %
        tpc_l{i}=(par_t+par_i)*(path{i}*lcost_l)';%计算测试的路径费用的拥堵费用部分
        mpc_l{i}=repmat(par_p.*(bench+(par_r-par_c).*ODflow_l{i}),p_size(i),1);%计算测试的路径费用的货币费用部分
        pcost_l{i}=tpc_l{i}(:)+mpc_l{i}+repmat(par_fix,p_size(i),1);%计算测试的路径费用
        dpf_l{i}=pflow{i}-pflow_l{i};
        dpc_l{i}=pcost{i}-pcost_l{i};
        al(i)=beta_l*dpc_l{i}'*dpc_l{i};
        be(i)=dpf_l{i}'*dpc_l{i};
    end
    
%下面判断步长    
    d_l=sum(al)/sum(be);
    if d_l<=2-delta
        d_rec(iter)=d_l;
        l_rec(iter)=l_l;%记录第iter次迭代的l_l
        beta_rec(iter)=beta_l;%记录第iter次迭代的步长beta_l
        beta=beta_l;
        pflow=pflow_l;
        pcost=pcost_l; 
        lcost=lcost_l;
        ODflow=ODflow_l;
        lflow=lflow_l;
        
        break;
    end
end
time(iter)=toc;
end
time_sum=sum(time);
    e=1;

    
%以下统计OD对的trips
trips=0;
for i=1:528
    o=OD(i,1);
    d=OD(i,2);    
    trips(o,d)=OD(i,3);
end

    
%以下统计OD对的路径数
sz1=0;sz2=0;sz3=0;sz4=0;sz5=0;sz6=0;sz7=0;sz_other=0;
for i=1:528
    o=OD(i,1);
    d=OD(i,2);    
    sz=size(H{i},1)/5;
    od_p_num(o,d)=sz;
    if sz==1
        sz1=sz1+1;
    elseif sz==2
        sz2=sz2+1;
    elseif sz==3
        sz3=sz3+1;
    elseif sz==4
        sz4=sz4+1;
    elseif sz==5
        sz5=sz5+1;
    elseif sz==6
        sz6=sz6+1;
    elseif sz==7
        sz7=sz7+1;
    else
        sz_other=sz_other+1;
    end
end
    sum(sum(od_p_num,1),2)
    
%以下生成初始约束矩阵
% for i=1:528
%         H{i}=eye(5*p_num);
%         p_beq{i}=zeros(2*p_num,1); 
%         beq{i}=[OD_demand(i);p_beq{i}];
%         R1(i,1:p_num) = {r1};
%         R2(i,1:p_num) = {r2};
%         A1{i} = blkdiag(R1{i,:});%生成流量约束1的系数矩阵（p_size维分块对角矩阵，每个分块等于[0,1,0,-1,0]）
%         A2{i} = blkdiag(R2{i,:});%生成流量约束2的系数矩阵（p_size维分块对角矩阵，每个分块等于[0,0,1,0,-0.5]）
%         d_Aeq{i}=ones(1,5*p_num);%需求约束的系数
%         Aeq{i}=[d_Aeq{i};A1{i};A2{i};];%等式约束的系数
%         lb{i}=zeros(5*p_num,1);%路径流量的下界
% end

for i=1:w
    path{i}=route_link((i-1)*p_num+1:i*p_num,:);
end

for i=1:w
    pflow{i}(1:5*p_num,1)=0;
end

%下面把OD流放到每个OD对的第一条路径
for i=1:w
    pflow{i}(1,1)=OD_demand(i);
end



%下面更新OD流量
for i=1:w
    B=reshape(pflow{i},5,[]);
    ODflow{i}=sum(B,2);    
end

%下面更新初始pcost费用
for i=1:w
    pcost{i}(1:5,1)=(par_t+par_i).*(path{i}(1,:)*link_cost)+par_p.*(bench+(par_r-par_c).*ODflow{i})+par_fix;
end

%计算每条路径上的流量
size_h=length(pflow);
size_l=link_number;
car_flow=zeros(size_h,1);
share_flow=zeros(size_h,1);
rider_flow=zeros(size_h,1);
car_flow1=zeros(size_h,1);
share_flow1=zeros(size_h,1);
rider_flow1=zeros(size_h,1);
for i=1:size_h
    car_flow(i)=pflow_l{i}(1)+pflow_l{i}(2)+pflow_l{i}(3);
    share_flow(i)=pflow_l{i}(2)+pflow_l{i}(3);
    rider_flow(i)=pflow_l{i}(4)+pflow_l{i}(5);
end
route_link_flow_1=zeros(size_h,size_l);
route_link_flow_2=zeros(size_h,size_l);
route_link_flow_3=zeros(size_h,size_l);
for i=1:size_h
    route_link_flow_1(i,:)=route_link(i,:)*car_flow(i);
    route_link_flow_2(i,:)=route_link(i,:)*share_flow(i);
    route_link_flow_3(i,:)=route_link(i,:)*rider_flow(i);
end
for i=1:size_h
    car_flow1=sum(route_link_flow_1,1);
    share_flow1=sum(route_link_flow_2,1);
    rider_flow1=sum(route_link_flow_3,1);
end
if ~exist('car_flow.xlsx','file')
    delete("car_flow.xlsx")
    delete('rider_flow.xlsx')
    delete('share_flow.xlsx')
end 
xlswrite('car_flow.xlsx',car_flow);
xlswrite('share_flow.xlsx',share_flow);
xlswrite('rider_flow.xlsx',rider_flow);