
1.����˴�����Ķ��Ͱ�Ӱ���취�������̵����

2.��ȫ�������������Ӧ��ע��

3.�������ʵĵط���readme.txt��ע�ͳ�����

�����ɻ�㣺

videoProcessing.m

   minValue(frame < minValue) = frame(frame < minValue); % ������������
   maxValue(frame > maxValue) = frame(frame > maxValue); 

   vRowPosEnter(idx) = (j-1) + (-vImg(idx,j-1))./(vImg(idx,j)-vImg(idx,j-1))+vCols(1)-1; % idx����������

   idx = (frame1 >= shadowValue) & (frame2 <  shadowValue) & isnan(shadowEnter); % ������ͼ�����ֵͼ��Ƚϣ��ڲ����������
   shadowEnter(idx) = (i-1) + (shadowValue(idx)-frame1(idx))./(frame2(idx)-frame1(idx)); % idx��ͼ�����У���

shadowScan.m
   
   X = [0 dX dX 0; 0 0 dY dY; 0 0 0 0];
   hPlane = fitPlane(X(1,:),X(2,:),X(3,:))'; %  [a b c d]';
   X = Rc_h'*(Rc_v*X + repmat(Tc_v-Tc_h,1,size(X,2)));
   vPlane = fitPlane(X(1,:),X(2,:),X(3,:))'; %  [a b c d]'; % ˮƽ��ʹ�ֱ����趨����������ѧ����

   shadowPlaneEnter(i,4) = 0.5*shadowPlaneEnter(i,1:3)*(q_v+q_h);% ͶӰ��ƽ��ķ������ĵ���ά

   ~isnan(shadowEnter(idx(i))) % ���� ~isnan �� idx(i)����ز���������
