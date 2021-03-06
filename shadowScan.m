% 重置matlab环境
clear; clc; close all;

% 添加工具包
addpath(genpath('./util'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 1: Define the data capture and reconstruction parameters.                                        
% Note: Uncomment a block to reconstruct the corresponding sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display progress.
disp('[Scanning with Shadows]');
disp(' ');  
disp('Loading object and reconstruction parameters...');

% Select which reconstruction script to use.
addpath('./data/man/');      man_v1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2: Video Processing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display progress.
disp('Performing video processing...');

% Modify extensions for low-resolution imagery (if enabled).
if useLowRes
   seqName = [seqName,'-lr'];
   calName = [calName,'-lr'];
end

% Get two points on the "middle" line (from user input).
disp('   + finding the dividing line between reference planes...');
disp('     (click on two points on the dividing line)');
frame = imread(['./data/',objName,'/',seqName,'/000001.jpg']);
figure(1); clf; set(gcf,'Name','Dividing Line Selection');
imagesc(frame); axis image;
title('Dividing Line Selection');
xlabel('click on two points on the dividing line');
middlePoints = ginput(2);

% Get bounding box for estimated "vertical" shadow boundary.
disp('   + define the reference area for the "vertical" plane...');
disp('     (click on top-left and bottom-right corners of reference area)');
figure(1); clf; set(gcf,'Name','Reference Area for "Vertical" Plane');
imagesc(frame); axis image;
title('Reference Area for "Vertical" Plane');
xlabel('click on the top-left and bottom-right corners of reference area');
x = ginput(2);
% 垂直面范围
vRows = round(min(x(:,2))):round(max(x(:,2))); % 如 vRows = 125 126 127
vCols = round(min(x(:,1))):round(max(x(:,1)));

% Get bounding box for estimated "horizontal" shadow boundary.
disp('   + define the reference area for the "horizontal" plane...');
disp('     (click on top-left and bottom-right corners of reference area)');
figure(1); clf; set(gcf,'Name','Reference Area for "Horizontal" Plane');
imagesc(frame); axis image;
title('Reference Area for "Horizontal" Plane');
title('Reference Area for "Horizontal" Plane');
xlabel('click on the top-left and bottom-right corners of reference area');
x = ginput(2);
% 水平面范围
hRows = round(min(x(:,2))):round(max(x(:,2)));
hCols = round(min(x(:,1))):round(max(x(:,1)));

% Perform necessary video processing.
disp('videoProcessing...');
% 见 videoProcessing.m
videoProcessing; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 3: Extrinsic Calibration.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display progress.
disp('Determining extrinsic calibration of reference plane(s)...');

% Load intrinsic calibration.
% Note: This assumes that you have already completed intrinsic calibration
%       for this sequence. See the assignment handout for more details.
load(['./data/',calName,'/Calib_Results'],'fc','cc','kc','alpha_c','nx','ny'); % 焦距、主点、畸变、偏斜、图像大小

% Obtain the extrinsic calibration for the "horizontal" plane.
% 需要手动选取各平面的四个基准点
% 获得水平平面的外部校正
disp('   + finding the extrinsic parameters of the "horizontal" plane....');
firstFrame = ['./data/',objName,'/',seqName,'/000001']; % ./data/man/v1-lr/000001
[Tc_h,Rc_h,H_h] = ... % 输出相机坐标系下水平面的平移向量、旋转矩阵、单应矩阵
    computeExtrinsic(fc,cc,kc,alpha_c,firstFrame,'jpg',nX,nY,dX,dY); %  dX,dY 标定板大小 558.8mm×303.2125mm
% [Tc_h,Rc_h,H_h]
% 输入焦距、主点、畸变、偏斜、图像名称、格式、单位像素、大小（毫米）
pause(1.5);

% Obtain the extrinsic calibration for the "vertical" plane.
% 获得垂直平面的外部校正
disp('   + finding the extrinsic parameters of the "vertical" plane....');
firstFrame = ['./data/',objName,'/',seqName,'/000001'];
[Tc_v,Rc_v,H_v] = ... % 输出相机坐标系下垂直面的平移向量、旋转矩阵、单应矩阵
   computeExtrinsic(fc,cc,kc,alpha_c,firstFrame,'jpg',nX,nY,dX,dY);
% [Tc_v,Rc_v,H_v]
pause(1.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 4: Reconstruction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display progress.
disp('Reconstructing 3D points using intersection with shadow plane(s)...');

% Estimate parameters of reference planes (using least-squares).
% 参考平面的估计参数(使用最小二乘法)
X = [0 dX dX 0; 0 0 dY dY; 0 0 0 0];
% 0  558.8000  558.8000         0
% 0         0  303.2125  303.2125
% 0         0         0         0
hPlane = fitPlane(X(1,:),X(2,:),X(3,:))'; %  [a b c d]';
%     0
%     0
%     1
%     0
X = Rc_h'*(Rc_v*X + repmat(Tc_v-Tc_h,1,size(X,2)));
%    1.2226  560.0225  559.8749    1.0750
%  403.8838  403.8154  401.4784  401.5468
%   95.5983   95.8698  399.0733  398.8018
vPlane = fitPlane(X(1,:),X(2,:),X(3,:))'; %  [a b c d]';
%    0.0001
%    1.0000
%    0.0077
%  404.6088

% Calculate camera center (in "horizontal" reference coordinate system).
C = -Rc_h'*Tc_h; % 2.4.4节公式 相机光心空间点，所有射线都过该点，相当于公式里的q

% Display progress.
disp('   + recovering implicit representation of shadow planes...');

% Determine implicit representation for the shadow planes.确定阴影面的隐式表示
% Note: There are several approaches for this step. See the IJCV paper
%       for a complete description, or consult the assignment handout.
shadowPlaneEnter = zeros(length(recFrames),4); % "entering" shadow plane
shadowPlaneLeave = zeros(length(recFrames),4); % "leaving" shadow plane
for i = 1:length(recFrames) 
   % 空间定位计算的是各顶点对应的空间点，并确定91张图片中的首光平面和尾光平面的方程，为插值提供依据

   % Determine true position of the "vertical" shadow boundary (entering).
   % intersectLines(vLineEnter(i,:),middleLine) 求两直线交点
   % =(x,y) =
   %         329.0472
   %         226.7120
   % pixel2ray(x,fc,cc,kc,alpha_c) % 计算相机射线方程
   %　upperLine　上交像素点
   %　middleLine　中间交像素点
   n1_v = Rc_h'*pixel2ray(intersectLines(vLineEnter(i,:),middleLine),fc,cc,kc,alpha_c); %　经过中间交点的相机光线方向向量
   % 0.9786
   % 0.1794
   % -0.1003
   n2_v = Rc_h'*pixel2ray(intersectLines(vLineEnter(i,:),upperLine),fc,cc,kc,alpha_c); %　经过上交点的相机光线方向向量
   % -0.9958
   % 0.0917
   % -0.0067
   % p = intersectLineWithPlane(q,v,w)  p = q + t*v;
   % finds the point of intersection　of a line in parametric form 
   % p = q + λv 
   %　with a plane W defined in implicit form.
   p1_v = intersectLineWithPlane(C,n1_v,vPlane);  
   % 经过光心C、方向向量为n1_v的射线与vPlane平面的交点p1_v
   % 9.5854
   % 0.4039
   % -0.0586
   p2_v = intersectLineWithPlane(C,n2_v,vPlane);  
   % 经过光心C、方向向量为n2_v的射线与vPlane平面的交点p2_v
   % -1.8009
   % 0.0401
   % 0.0765
   
   % Determine true position of the "horizontal" shadow boundary (entering).
   n1_h = Rc_h'*pixel2ray(intersectLines(hLineEnter(i,:),middleLine),fc,cc,kc,alpha_c); %　经过中间交点的相机光线方向向量
   n2_h = Rc_h'*pixel2ray(intersectLines(hLineEnter(i,:),lowerLine),fc,cc,kc,alpha_c);  %　经过下交点的相机光线方向向量
   p1_h = intersectLineWithPlane(C,n1_h,hPlane); 
   % 经过光心C、方向向量为n1_h的射线与hPlane平面的交点p1_h
   p2_h = intersectLineWithPlane(C,n2_h,hPlane); 
   % 经过光心C、方向向量为n2_h的射线与hPlane平面的交点p2_h
   
   % Compute the "entering" shadow plane parameters.
   % 计算边缘的空间方程
   q_v = p1_v;  % 首边缘空间起点(中间交点) 
   v_v = (p2_v-p1_v)/norm(p2_v-p1_v); % 方向向量 （从中间交点指向上交点）
   q_h = p1_h;  % 首边缘空间起点(中间交点) 
   v_h = (p2_h-p1_h)/norm(p2_h-p1_h); % 方向向量 （从中间交点指向下交点）
   shadowPlaneEnter(i,1:3) = cross(v_v,v_h);
   % 0.2098   -0.5543    0.8054
   % 求投影光平面的法向量，前三维由叉积决定 
   % C = CROSS(A,B)  C = A x B.  A and B must be 3 element vectors.
   shadowPlaneEnter(i,1:3) = shadowPlaneEnter(i,1:3)/norm(shadowPlaneEnter(i,1:3)); % 归一化
   shadowPlaneEnter(i,4) = 0.5*shadowPlaneEnter(i,1:3)*(q_v+q_h); % 如果不相交就取中点，看文献
   % -6.5581
   
   % Determine true position of the "vertical" shadow boundary (leaving).
   n1_v = Rc_h'*pixel2ray(intersectLines(vLineLeave(i,:),middleLine),fc,cc,kc,alpha_c);
   n2_v = Rc_h'*pixel2ray(intersectLines(vLineLeave(i,:),upperLine),fc,cc,kc,alpha_c);
   p1_v = intersectLineWithPlane(C,n1_v,vPlane);
   p2_v = intersectLineWithPlane(C,n2_v,vPlane);
   
   % Determine true position of the "horizontal" shadow boundary (leaving).
   n1_h = Rc_h'*pixel2ray(intersectLines(hLineLeave(i,:),middleLine),fc,cc,kc,alpha_c);
   n2_h = Rc_h'*pixel2ray(intersectLines(hLineLeave(i,:),lowerLine),fc,cc,kc,alpha_c);
   p1_h = intersectLineWithPlane(C,n1_h,hPlane);
   p2_h = intersectLineWithPlane(C,n2_h,hPlane);
   
   % Compute the "leaving" shadow plane parameters.
   q_v = p1_v;
   v_v = (p2_v-p1_v)/norm(p2_v-p1_v);
   q_h = p1_h;
   v_h = (p2_h-p1_h)/norm(p2_h-p1_h);
   shadowPlaneLeave(i,1:3) = cross(v_v,v_h);
   shadowPlaneLeave(i,1:3) = shadowPlaneLeave(i,1:3)/norm(shadowPlaneLeave(i,1:3));
   shadowPlaneLeave(i,4) = 0.5*shadowPlaneLeave(i,1:3)*(q_v+q_h);
   
end

% Display progress.
disp('   + reconstructing 3D points...');

% Load first frame for assigning color.
% Note: You could use the maximum value image or some other source here.
frame = im2double(imread(['./data/',objName,'/',seqName,'/',num2str(allFrames(1),'%0.6d'),'.jpg']));

% Reconstruct 3D points using intersection with shadow plane(s).
% Note: If multiple shadow planes are used, then reduce errors
%       by averaging available estimates. This can also assist
%       in removing some outliers (if estimates do not agree).
% a = ~isnan(shadowEnter); 384X512 logical
% b = ~isnan(shadowLeave); 384X512 logical
% c = a & b; 384X512 logical
idx       = find(~isnan(shadowEnter) & ~isnan(shadowLeave)); % 先取逻辑值，再找索引 13289X1 double
% find（A）返回矩阵A中非零元素所在位置,从第一行开始，1-inf
% ~isnan(x)求得逻辑值
% 13289X1 double
idx       = idx(1:dSample:length(idx)); % 重采样，采样率1/5
[row,col] = ind2sub(size(shadowEnter),idx);% 把数组中元素索引值转换为该元素在数组中对应的下标，即扫描点
%   358   262
%   363   262
%   368   262
%   ...
npts      = length(idx);
vertices  = zeros(npts,3);
colors    = 0.65*ones(npts,3);
h = waitbar(0,'Reconstructing 3D points...');
for i = 1:npts % 被扫描的点有13289个，依次计算其空间位置
   
   % Obtain the camera ray for this pixel.
   n = Rc_h'*pixel2ray([col(i); row(i)],fc,cc,kc,alpha_c);
   %　经过物体表面的像素点（col(i),row(i)）的相机光线方向向量
   % 0.2499
   % 0.8134
   % -0.5253

   % Interpolate "entering" shadow plane parameters (using shadow time).
   % a = idx(i);
   % b = shadowEnter(idx(i));
   % c = ~isnan(shadowEnter(idx(i)));
   % a =
   %    170075
   % b =
   %    6.9688
   % c =
   %    1
   if ~isnan(shadowEnter(idx(i)))
      t  = shadowEnter(idx(i)); % 相当于第t幅图像
      t1 = floor(t); % 向下取整
      t2 = t1+1;
      if t2 <= length(recFrames)
         alpha = (t-t1)/(t2-t1); % 比例
         wEnter = (1-alpha)*shadowPlaneEnter(t1,:)+alpha*shadowPlaneEnter(t2,:); 
         % shadowPlaneEnter：length(recFrames)X4  91X4 
         % 按比例插值计算在第t1幅光平面和第t2幅光平面之间的wEnter光平面方程
         pEnter = intersectLineWithPlane(C,n,wEnter'); % 经过光心C、方向向量为n的射线与wEnter平面的交点pEnter
         vertices(i,:) = pEnter; % 物体表面像素点对应的空间点
      end
   end
   
   % Interpolate "leaving" shadow plane parameters (using shadow time).
   if ~isnan(shadowLeave(idx(i)))
      t  = shadowLeave(idx(i));
      t1 = floor(t);
      t2 = t1+1;
      if t2 <= length(recFrames)
         alpha = (t-t1)/(t2-t1);
         wLeave = (1-alpha)*shadowPlaneLeave(t1,:)+alpha*shadowPlaneLeave(t2,:);
         pLeave = intersectLineWithPlane(C,n,wLeave');
         vertices(i,:) = pLeave;
      end
   end

   % Average "entering" and "leaving" estimates (if both are available).
   % Note: If points do not agree, set coordinate to infinity.
   %       This will ensure that it is clipped by the bounding volume.
   % 如果点不一致，设坐标为无穷。这将确保它被边界卷所限。
   if ~isnan(shadowEnter(idx(i))) && ~isnan(shadowLeave(idx(i)))
      if norm(pEnter-pLeave) <= distReject % 首尾距离
         vertices(i,:) = 0.5*(pEnter + pLeave);
      else
         vertices(i,:) = Inf*[1 1 1];
      end
   end
   
   % Assign color per vertex (using source image).
   colors(i,:) = frame(row(i),col(i),:); % 着色
   
   % Update progress bar.
   waitbar(i/npts);
      
end
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 5: Post-processing and Visualization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display progress.
disp('Display reconstruction results and exporting VRML file...');

% Display the 3D configuraton using the extrinsic calibration.
% Note: Some of this is hard-coded for the configuration included
%       with the assignment. This could be generalized.
disp('   + displaying reconstruction results...');
figure(1); set(gcf,'Name','Reconstruction Results'); clf;
X = [0 dX dX 0; 0 0 dY dY; 0 0 0 0];
patch(X(1,:),X(2,:),X(3,:),0.65*[1 1 1],'FaceAlpha',0.5);
hold on;
   % plots the "horizontal" plane
   plot3(X(1,:),X(2,:),X(3,:),'g.','MarkerSize',15); 
hold off;
X = Rc_h'*(Rc_v*X + repmat(Tc_v-Tc_h,1,size(X,2)));
patch(X(1,:),X(2,:),X(3,:),0.65*[1 1 1],'FaceAlpha',0.5);
hold on;
   % plots the "vertical" plane
   plot3(X(1,:),X(2,:),X(3,:),'r.','MarkerSize',15);
hold off;
C = -Rc_h'*Tc_h;
hold on;
   % plots the camera center
   plot3(C(1),C(2),C(3),'b.','MarkerSize',20);
hold off;
axis equal tight; grid on; view(3); box on;
axis([-200 800 -1300 500 -100 1000]);
xlabel('x'); ylabel('y'); zlabel('z');

% Clip the recovered 3D point cloud using the bounding volume.
clip = find( (vertices(:,1) >= clipRangeX(1) & vertices(:,1) <= clipRangeX(2)) & ...
             (vertices(:,2) >= clipRangeY(1) & vertices(:,2) <= clipRangeY(2)) & ...
             (vertices(:,3) >= clipRangeZ(1) & vertices(:,3) <= clipRangeZ(2)) );

% Display the recovered 3D point cloud (with per-vertex color).
% Note: Convert to indexed color map for use with FSCATTER3.
figure(1); set(gcf,'Name','Reconstruction Results');
C = reshape(colors,[size(colors,1) 1 size(colors,2)]);
[C,cmap] = rgb2ind(C,256);
hold on;
   %plot3(vertices(:,1),vertices(:,2),vertices(:,3),'b.','MarkerSize',5);
   fscatter3(vertices(:,1),vertices(:,2),vertices(:,3),double(C),cmap);
hold off;

% Export colored point cloud as a VRML file.
% Note: Interchange x and y coordinates for j3DPGP.
disp('   + exporting VRML file...');
writeVRML(['./models/',objName,'_',seqName,'.wrl'],...
   vertices(clip,[2 1 3]),...
   colors(clip,:));
disp(' ');
