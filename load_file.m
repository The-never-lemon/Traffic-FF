function [name,p_name]=load_file()
[Fnameh,Pnameh]=uigetfile('*.mat');%打开文件浏览器选择.mat文件
 file_name = strcat(Pnameh,Fnameh);%获取包含地址的文件名
 name=file_name;
 p_name=Pnameh;
end
