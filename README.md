# Traffic-FF: A traffic assignment software considering multiple ridesharing services
## 1. Instructions
Traffic-FF is an application used to predict traffic flow. It's developed based on **MATLAB** and **MATLAB APP Designer**.
## 2. Installation
* __Necessary toolboxes__  
`Parallel Computing Toolbox`  
`Optimization Toolbox`  
`MATLAB Compiler`(if you want to package the app)
* __Download the code__  
You can download the code from this repository including `dijkstra.m` `kShortestPath.m` `load_file.m` `parallel_column_sf_1030.m` `sfKShortest.m`.
After putting them in the same path,you just need to run `parallel_column_sf_1030.m`.   
* __Download the app__  
Apart from downloading the code to run, you can also use the packaged app.The whole steps are shown in the [oi-traffic-FF.pdf](./oi-Traffic-FF.pdf).
## 3. Test data
Test data is evolved from the topology of the city of Sioux-Falls in the USA, which you can find in [column1_0323_2.5OD.mat](./column1_0323_2.5OD.mat).The Sioux-Falls network has 24 nodes, 76 links, and 528 OD pairs.The mat contains the links and OD pairs of the city
The detailed parameter values of the links and OD pairs can be found in [Bar-Gera](https://github.com/Gerald-Development/Barista-Gerald).
## 4. Contact us 
If you find problems or bugs, you can leave your comments and send us emails. Here are the emails.
* majie@seu.edu.cn
* jrzhu@seu.edu.cn
