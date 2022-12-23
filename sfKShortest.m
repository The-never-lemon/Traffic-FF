clear;
clc;
load('CGsf2.mat');%mat�ں������������-OD��link��route_link��OD_route
%������
link_number=size(link,1);%�õ�·����
OD_number=size(OD,1);%�õ�OD��
OD_demand=zeros(OD_number,1,'double');%OD������
OD_route_number=zeros(OD_number,1,'double');%ÿ��OD��Ӧ��·����
link_volume=zeros(link_number,1,'double');%·������
link_volume_pre=zeros(link_number,1,'double');%ǰһ�ε�·������
link_cost=zeros(link_number,1,'double');%·�η���
d_link_cost=zeros(link_number,1,'double');%·�η��õĵ���
link_capacity=zeros(link_number,1,'double');%·��ͨ������
link_t0=zeros(link_number,1,'double');%·������������ʱ��
link_b0=zeros(link_number,1,'double');%·������Լ����Ӧ���������ճ���
link_r0=zeros(link_number,1,'double');%��ǿ�������ճ��ӷ��еĳͷ�����
link_origin=zeros(link_number,1,'double');%·�ε����
link_destination=zeros(link_number,1,'double');%·�ε��յ�
link_zero=zeros(link_number,1,'double');%���ڼ����������ճ���
link_Ca_cost=zeros(link_number,1,'double');%·���ڴﵽͨ������ʱ�ķ���
route_link_1=zeros(1,link_number,'double');%·���ڴﵽͨ������ʱ�ķ���
x_change=zeros(link_number,1,'double');%·������ÿ�εı仯��
min_route=zeros(OD_number,1,'double');%�洢ÿ�ε���������ĸ���OD�Ե����·��
min_route_cost=zeros(OD_number,1,'double');%�洢ÿ�ε���������ĸ���OD�Ե����·��

OD_demand=OD(:,3);%OD������
link_t0=link(:,3);%·������������ʱ��
link_cost=link_t0;%·�η���
link_capacity=link(:,4)*1000;%·��ͨ������
link_origin=link(:,1);%·������������ʱ��;%·�ε����
link_destination=link(:,2);%·������������ʱ��;%·�ε��յ�

%�����·
DG=sparse(link_origin',link_destination',link_cost);
for i=1:length(DG)
    for j=1:length(DG)
        if DG(i,j)==0
            DG(i,j)=Inf;
        end
    end
end
route_num=1;%·���Ǻ�
for i=1:length(OD_demand)
[routeCell,Dist]=kShortestPath(DG,OD(i,1),OD(i,2),1);%���·�㷨�����������[k�����·�����Խ��Ϊ��ǣ�cell��,����]
    i_path_num=1;
    for ci =1:length(routeCell)
        route=cell2mat(routeCell(ci));%�ѵ�ci��·���Ľ���ʾ��mat
        for j=1:size(route,2)-1
            [link_ID,l]=find(link(:,5)==route(j)*100+route(j+1));%����·��Ϊ��ǵ����·ת�����Թ�������Ϊ���
            route_link(route_num,link_ID)=1;%���ɹ�������
        end        
        OD_route(i,route_num)=1;
        route_num=route_num+1;
        
    end
    %OD_route_number(i)=route_num;
end

