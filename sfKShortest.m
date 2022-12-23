clear;
clc;
load('CGsf2.mat');%mat内含各种输入矩阵-OD、link、route_link、OD_route
%主程序
link_number=size(link,1);%得到路段数
OD_number=size(OD,1);%得到OD数
OD_demand=zeros(OD_number,1,'double');%OD需求量
OD_route_number=zeros(OD_number,1,'double');%每个OD对应的路径数
link_volume=zeros(link_number,1,'double');%路段流量
link_volume_pre=zeros(link_number,1,'double');%前一次的路段流量
link_cost=zeros(link_number,1,'double');%路段费用
d_link_cost=zeros(link_number,1,'double');%路段费用的导数
link_capacity=zeros(link_number,1,'double');%路段通行能力
link_t0=zeros(link_number,1,'double');%路段自由流旅行时间
link_b0=zeros(link_number,1,'double');%路段容量约束对应的拉格朗日乘子
link_r0=zeros(link_number,1,'double');%增强拉格朗日乘子法中的惩罚因子
link_origin=zeros(link_number,1,'double');%路段的起点
link_destination=zeros(link_number,1,'double');%路段的终点
link_zero=zeros(link_number,1,'double');%用于计算拉格朗日乘子
link_Ca_cost=zeros(link_number,1,'double');%路段在达到通行能力时的费用
route_link_1=zeros(1,link_number,'double');%路段在达到通行能力时的费用
x_change=zeros(link_number,1,'double');%路段流量每次的变化量
min_route=zeros(OD_number,1,'double');%存储每次迭代计算出的各个OD对的最短路径
min_route_cost=zeros(OD_number,1,'double');%存储每次迭代计算出的各个OD对的最短路径

OD_demand=OD(:,3);%OD需求量
link_t0=link(:,3);%路段自由流旅行时间
link_cost=link_t0;%路段费用
link_capacity=link(:,4)*1000;%路段通行能力
link_origin=link(:,1);%路段自由流旅行时间;%路段的起点
link_destination=link(:,2);%路段自由流旅行时间;%路段的终点

%求最短路
DG=sparse(link_origin',link_destination',link_cost);
for i=1:length(DG)
    for j=1:length(DG)
        if DG(i,j)==0
            DG(i,j)=Inf;
        end
    end
end
route_num=1;%路径记号
for i=1:length(OD_demand)
[routeCell,Dist]=kShortestPath(DG,OD(i,1),OD(i,2),1);%最短路算法，求出来的是[k条最短路径（以结点为标记，cell）,费用]
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

