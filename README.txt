
1.完成了代码的阅读和暗影铸造法整体流程的理解

2.对全部代码进行了相应的注释

3.对有疑问的地方在readme.txt中注释出来了

程序疑惑点：

videoProcessing.m

   minValue(frame < minValue) = frame(frame < minValue); % 具体的运算过程
   maxValue(frame > maxValue) = frame(frame > maxValue); 

   vRowPosEnter(idx) = (j-1) + (-vImg(idx,j-1))./(vImg(idx,j)-vImg(idx,j-1))+vCols(1)-1; % idx是行数？？

   idx = (frame1 >= shadowValue) & (frame2 <  shadowValue) & isnan(shadowEnter); % 相邻两图像和阈值图像比较，内部的运算过程
   shadowEnter(idx) = (i-1) + (shadowValue(idx)-frame1(idx))./(frame2(idx)-frame1(idx)); % idx是图像序列？？

shadowScan.m
   
   X = [0 dX dX 0; 0 0 dY dY; 0 0 0 0];
   hPlane = fitPlane(X(1,:),X(2,:),X(3,:))'; %  [a b c d]';
   X = Rc_h'*(Rc_v*X + repmat(Tc_v-Tc_h,1,size(X,2)));
   vPlane = fitPlane(X(1,:),X(2,:),X(3,:))'; %  [a b c d]'; % 水平面和垂直面的设定参数？？数学基础

   shadowPlaneEnter(i,4) = 0.5*shadowPlaneEnter(i,1:3)*(q_v+q_h);% 投影光平面的法向量的第四维

   ~isnan(shadowEnter(idx(i))) % 包含 ~isnan 和 idx(i)等相关参数的运算
