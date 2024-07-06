clc
close all
clear all
syms alpha theta phi sigma real
 set(gcf, 'color', 'white');
set(gca,'visible','off')
% kinematic parameters % dimensions in mm
% spherical joint base home pose
base_s=[0;0;0];
%  p_os=[0;0;72];
% b_s1_o=[0;129;0];b_s2_o=[0;129;0];b_s3_o=[-111.7;-64.5;0];b_s4_o=[-111.7;-64.5;0];
% b_s5_o=[111.7;-64.5;0];b_s6_o=[111.7;-64.5;0];
% p_s1_s=[42.4;24.5;51];p_s6_s=[42.4;24.5;51];p_s2_s=[-42.4;24.5;51];p_s3_s=[-42.4;24.5;51];
% p_s4_s=[0;-49;51];p_s5_s=[0;-49;51];

% given below are set of parameters in meters only for computation of
% polytope available wrench set
 p_os=0.001*[0;0;72];
% p_os=0.001*[0;0;0];
b_s1_o=0.001*[0;129;0];b_s2_o=0.001*[0;129;0];b_s3_o=0.001*[-111.7;-64.5;0];b_s4_o=0.001*[-111.7;-64.5;0];
b_s5_o=0.001*[111.7;-64.5;0];b_s6_o=0.001*[111.7;-64.5;0];
p_s1_s=0.001*[42.4;24.5;51];p_s6_s=0.001*[42.4;24.5;51];p_s2_s=0.001*[-42.4;24.5;51];p_s3_s=0.001*[-42.4;24.5;51];
p_s4_s=0.001*[0;-49;51];p_s5_s=0.001*[0;-49;51];
shadow=0;%1=shadow, 0=none
% shoulder joint limit phi [-180,180], theta [0,140] sigma [-90,90]
%   set(gcf, 'color', 'white','units','pixels','position',[0 0 1920 1080]);

% vidObj=VideoWriter('4spsshoulder3.mp4');
% vidObj.FrameRate = 24;
% open(vidObj);
t_ini=1.5;
% t_fin=2.6;
t_fin=6.6;
R=2;


% phi=   -2.6175;
% theta=-0.5236;
% sigma=    -5.7589;
%    sigma=deg2rad(47);
% for t1=t_ini
%  for t1=t_ini:0.1:t_fin
     syms t real
%     xt=R*cos(0.5*t) 
%     yt=R*cos(0.5*t)
%     zt= zeros(1,numel(xt))+123;
%     fplot(xt,yt)
%     fplot3(xt,yt,zt)
%  theta=R*cos(0.5*t1)
%  phi=R*sin(0.5*t1)

% phi=0;  % phi theta psi
%  theta=pi/2;
%   sigma=0;
phi=deg2rad(15);  % phi theta psi
 theta=deg2rad(0);
  sigma=deg2rad(15);
%  phi=deg2rad(40);  % phi theta psi
%  theta=deg2rad(-50);
%   sigma=deg2rad(40);
 
% hold on
%  phi=0;  
% sigma=0;
% for phi=0:pi
%  for theta=pi/2:0.1:deg2rad(100)
%      for phi=deg2rad(80)
%          for sigma=0:0.1:pi/2
%kinematics of 3 dof cdpm with passive spherical joint shoulder
%   R_os=Rz(phi)*Ry(theta)*Rz(sigma-phi); %tilt and torsion, modified ZYZ
 R_os=[cos(phi)*cos(theta) cos(phi)*sin(theta)*sin(sigma)+sin(phi)*cos(sigma) -cos(phi)*sin(theta)*cos(sigma)+sin(phi)*sin(sigma)
     -sin(phi)*cos(theta) -sin(phi)*sin(theta)*sin(sigma)+cos(phi)*cos(sigma) sin(phi)*sin(theta)*cos(sigma)+cos(phi)*sin(sigma)
     sin(theta)  -cos(theta)*sin(sigma)  cos(theta)*cos(sigma)];
% %  R_os=Rx(phi)*Ry(theta)*Rz(sigma);
%  R_os=Rz(phi)*Ry(theta)*Rz(sigma); %ZYZ Euler
 R_os=Rz(phi)*Ry(theta)*Rx(sigma-phi);
T_os=[R_os p_os
    zeros(1,3) 1];

%shoulder
l_s1= norm(mtimes(R_os,p_s1_s)-b_s1_o)
l_s2= norm(mtimes(R_os,p_s2_s)-b_s2_o)
l_s3= norm(mtimes(R_os,p_s3_s)-b_s3_o)
l_s4= norm(mtimes(R_os,p_s4_s)-b_s4_o)
l_s5= norm(mtimes(R_os,p_s5_s)-b_s5_o)
l_s6= norm(mtimes(R_os,p_s6_s)-b_s6_o)

%base triangle 
xs_b= [b_s1_o(1) b_s3_o(1) b_s5_o(1)];% here triangle is formed with three unique vertices, other three points are same here otherwise need to plot polygon
ys_b= [b_s1_o(2) b_s3_o(2) b_s5_o(2)];
zs_b= [b_s1_o(3) b_s3_o(3) b_s5_o(3)];
fill3(xs_b,ys_b,zs_b,'r','LineStyle','none')
 material shiny
hold on 
%transformed coordinates
p_s1_o=R_os*p_s1_s+p_os;
p_s2_o=R_os*p_s2_s+p_os;
p_s3_o=R_os*p_s3_s+p_os;
p_s4_o=R_os*p_s4_s+p_os;
p_s5_o=R_os*p_s5_s+p_os;
p_s6_o=R_os*p_s6_s+p_os;
% xs_s= [p_s1_s(1) p_s3_s(1) p_s5_s(1)];% here triangle is formed with three unique vertices, other three points are same here otherwise need to plot polygon
% ys_s= [p_s1_s(2) p_s3_s(2) p_s5_s(2)];
% zs_s= [p_s1_s(3) p_s3_s(3) p_s5_s(3)];
xs_s= [p_s1_o(1) p_s3_o(1) p_s5_o(1)];% here triangle is formed with three unique vertices, other three points are same here otherwise need to plot polygon
ys_s= [p_s1_o(2) p_s3_o(2) p_s5_o(2)];
zs_s= [p_s1_o(3) p_s3_o(3) p_s5_o(3)];
%  axes('nextplot', 'add');
fill3(xs_s,ys_s,zs_s,'r','LineStyle','none')
% material shiny
%  axes('nextplot', 'add');
drawplanethreepoint(p_s1_o,p_s3_o,p_s5_o,'red','none',0) % first facecolor, second edgecolor
hold on
%133
%355
%151
[normal1,d1]=drawplanethreepoint(b_s1_o,b_s3_o,p_s3_o,'none','black',0) % first facecolor, second edgecolor
hold on
[normal2,d2]=drawplanethreepoint(b_s3_o,b_s5_o,p_s5_o,'none','black',0) % first facecolor, second edgecolor
hold on
[normal3,d3]=drawplanethreepoint(b_s1_o,b_s5_o,p_s1_o,'none','black',0) % first facecolor, second edgecolor
hold on
P=threeplaneintersectionpoint(normal1, normal2,normal3,d1,d2,d3)
plot3(P(1),P(2),P(3),'o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')
plot3([P(1) b_s1_o(1)]',[P(2) b_s1_o(2)]',[P(3) b_s1_o(3)]','Color' , 'k', 'LineWidth', 2); %ax1
hold on
plot3([P(1) b_s3_o(1)]',[P(2) b_s3_o(2)]',[P(3) b_s3_o(3)]','Color' , 'k', 'LineWidth', 2); %ax1
hold on
plot3([P(1) b_s5_o(1)]',[P(2) b_s5_o(2)]',[P(3) b_s5_o(3)]','Color' , 'k', 'LineWidth', 2); %ax1
hold on

hold on
%  axes('nextplot', 'add');
center_plat= [sum(xs_s)/3 sum(ys_s)/3 sum(zs_s)/3];
jointfc='k';
linkec='k';
hold on
% [x0,y0,z0]=box2P(1.7,[base_s(1),base_s(2),base_s(3)],[p_os(1),p_os(2),p_os(3)]);
% surf(x0,y0,z0,'FaceColor','r','Facealpha',0.6, ...
% 'EdgeColor',linkec,'AmbientStrength',0.6)
% patch('Vertices',[x0(1,:)' y0(1,:)' z0(1,:)'], ...
% 'Faces',[1:5],'FaceColor',jointfc, ...
% 'EdgeColor','none','AmbientStrength',0.6)
% patch('Vertices',[x0(2,:)' y0(2,:)' z0(2,:)'], ...
% 'Faces',[1:5],'FaceColor',jointfc, ...
% 'EdgeColor','none','AmbientStrength',0.6)
% hold on
% [x0,y0,z0]=box2P(1.7,[center_plat(1),center_plat(2),center_plat(3)],[p_os(1),p_os(2),p_os(3)]);
% surf(x0,y0,z0,'FaceColor','r','Facealpha',0.6, ...
% 'EdgeColor',linkec,'AmbientStrength',0.6)
% patch('Vertices',[x0(1,:)' y0(1,:)' z0(1,:)'], ...
% 'Faces',[1:5],'FaceColor',jointfc, ...
% 'EdgeColor','none','AmbientStrength',0.6)
% patch('Vertices',[x0(2,:)' y0(2,:)' z0(2,:)'], ...
% 'Faces',[1:5],'FaceColor',jointfc, ...
% 'EdgeColor','none','AmbientStrength',0.6)
%% if scaled, scaled version that is in meters
[x0,y0,z0]=box2P(0.003,[base_s(1),base_s(2),base_s(3)],[p_os(1),p_os(2),p_os(3)]);
surf(x0,y0,z0,'FaceColor','r','Facealpha',0.6, ...
'EdgeColor',linkec,'AmbientStrength',0.6)
patch('Vertices',[x0(1,:)' y0(1,:)' z0(1,:)'], ...
'Faces',[1:5],'FaceColor',jointfc, ...
'EdgeColor','none','AmbientStrength',0.6)
patch('Vertices',[x0(2,:)' y0(2,:)' z0(2,:)'], ...
'Faces',[1:5],'FaceColor',jointfc, ...
'EdgeColor','none','AmbientStrength',0.6)
hold on
[x0,y0,z0]=box2P(0.003,[center_plat(1),center_plat(2),center_plat(3)],[p_os(1),p_os(2),p_os(3)]);
surf(x0,y0,z0,'FaceColor','r','Facealpha',0.6, ...
'EdgeColor',linkec,'AmbientStrength',0.6)
patch('Vertices',[x0(1,:)' y0(1,:)' z0(1,:)'], ...
'Faces',[1:5],'FaceColor',jointfc, ...
'EdgeColor','none','AmbientStrength',0.6)
patch('Vertices',[x0(2,:)' y0(2,:)' z0(2,:)'], ...
'Faces',[1:5],'FaceColor',jointfc, ...
'EdgeColor','none','AmbientStrength',0.6)
%% plot joint to platform lines just for visualization
lh1 =plot3([p_os(1) b_s1_o(1)]',[p_os(2) b_s1_o(2)]',[p_os(3) b_s1_o(3)]','Color' , 'r', 'LineWidth', 1, 'LineStyle',':'); %ax2
lh1.Color = [lh1.Color 0.5];
hold on
lh2 =plot3([p_os(1) b_s3_o(1)]',[p_os(2) b_s3_o(2)]',[p_os(3) b_s3_o(3)]','Color' , 'r', 'LineWidth', 1, 'LineStyle',':'); %ax2
lh2.Color = [lh2.Color 0.5];
hold on
lh3 =plot3([p_os(1) b_s5_o(1)]',[p_os(2) b_s5_o(2)]',[p_os(3) b_s5_o(3)]','Color' , 'r', 'LineWidth', 1, 'LineStyle',':'); %ax2
lh3.Color = [lh3.Color 0.5];
hold on
% top tetrahedron 
lh1 =plot3([p_os(1) p_s1_o(1)]',[p_os(2) p_s1_o(2)]',[p_os(3) p_s1_o(3)]','Color' , 'r', 'LineWidth', 1, 'LineStyle',':'); %ax2
lh1.Color = [lh1.Color 0.5];
hold on
lh2 =plot3([p_os(1) p_s3_o(1)]',[p_os(2) p_s3_o(2)]',[p_os(3) p_s3_o(3)]','Color' , 'r', 'LineWidth', 1, 'LineStyle',':'); %ax2
lh2.Color = [lh2.Color 0.5];
hold on
lh3 =plot3([p_os(1) p_s5_o(1)]',[p_os(2) p_s5_o(2)]',[p_os(3) p_s5_o(3)]','Color' , 'r', 'LineWidth', 1, 'LineStyle',':'); %ax2
lh3.Color = [lh3.Color 0.5];
hold on
%% plot joints
facecolor='none';
    edgecolor='k';
   jointfc='k';
%  jointfc= [0.6350 0.0780 0.1840];
%  jointec= [0.6350 0.0780 0.1840];
    jointec='k';
linkfc='r';
 linkfc2='r';
linkec='k';
k1_B=p_os-base_s/norm(p_os-base_s);
% sphere_joint_axis(10,k1_B,3.5,eye(3),p_os',jointfc,jointec)
sphere_joint_axis(10,k1_B,0.006,eye(3),p_os',jointfc,jointec)
hold on
light('style','local','position',[100 100 150],'color',[1 1 1])
%  light('style','local','position',[0 0 3],'color',[1 1 1])
 box on
%       daspect([1 1 1])
%     axis([xmin xmax ymin ymax zmin zmax])
%      axis vis3d
%% shadow
%Camera normal vector
    V1=xs_s;
    V2=ys_s;
    normal=cross(V1,V2)/norm(cross(V1,V2));
%     camtarget([0 0 .7])
%             camzoom(1.5)
%Plotting shadow
switch shadow
  case 1
shadowcolor=.95*[1 1 1];
patch(xa,ya,shadowcolor,'edgecolor','none')
otherwise
end
    hold on       
% %% lines
% plotlineparallel(b_s1_o,b_s3_o,p_s1_o)
% hold on
 %% points
 plot3(b_s1_o(1),b_s1_o(2),b_s1_o(3),'o','MarkerFaceColor','k')
 hold on
 plot3(b_s3_o(1),b_s3_o(2),b_s3_o(3),'o','MarkerFaceColor','k')
 hold on
 plot3(b_s5_o(1),b_s5_o(2),b_s5_o(3),'o','MarkerFaceColor','k')
 hold on
 
 plot3(p_s1_o(1),p_s1_o(2),p_s1_o(3),'o','MarkerFaceColor','k')
 hold on
 plot3(p_s3_o(1),p_s3_o(2),p_s3_o(3),'o','MarkerFaceColor','k')
 hold on
 plot3(p_s5_o(1),p_s5_o(2),p_s5_o(3),'o','MarkerFaceColor','k')
 hold on
% %% plot plane
% drawplanethreepoint(p_s1_o,p_s3_o,p_s5_o)
% hold on
%% adding labels
% text(b_s1_o(1)-0.01,b_s1_o(2)-0.01,b_s1_o(3),'bs1','fontsize',14,'Interpreter', 'latex');
% hold on
% text(b_s3_o(1)-0.01,b_s3_o(2)-0.01,b_s3_o(3),'bs3','fontsize',14,'Interpreter', 'latex');
% hold on
% text(b_s5_o(1)-0.01,b_s5_o(2)-0.01,b_s5_o(3),'bs5','fontsize',14,'Interpreter', 'latex');
% hold on
% text(p_s1_o(1)-0.01,p_s1_o(2)-0.01,p_s1_o(3),'ps1','fontsize',14,'Interpreter', 'latex');
% hold on
% text(p_s3_o(1)-0.01,p_s3_o(2)-0.01,p_s3_o(3),'ps3','fontsize',14,'Interpreter', 'latex');
% hold on
% text(p_s5_o(1)-0.01,p_s5_o(2)-0.01,p_s5_o(3),'ps5','fontsize',14,'Interpreter', 'latex');
% hold on

 %% plot cables
 jointfc='r';
 jointec='k';
 linkfc='r';
   prismatic_joint_axis(0.0015,[b_s1_o(1),b_s1_o(2),b_s1_o(3)],[p_s1_o(1),p_s1_o(2),p_s1_o(3)],jointfc,jointec,linkfc);
 hold on 
  prismatic_joint_axis(0.0015,[b_s2_o(1),b_s2_o(2),b_s2_o(3)],[p_s2_o(1),p_s2_o(2),p_s2_o(3)],jointfc,jointec,linkfc);
 hold on
  prismatic_joint_axis(0.0015,[b_s3_o(1),b_s3_o(2),b_s3_o(3)],[p_s3_o(1),p_s3_o(2),p_s3_o(3)],jointfc,jointec,linkfc);
 hold on
  prismatic_joint_axis(0.0015,[b_s4_o(1),b_s4_o(2),b_s4_o(3)],[p_s4_o(1),p_s4_o(2),p_s4_o(3)],jointfc,jointec,linkfc);
 hold on
   prismatic_joint_axis(0.0015,[b_s5_o(1),b_s5_o(2),b_s5_o(3)],[p_s5_o(1),p_s5_o(2),p_s5_o(3)],jointfc,jointec,linkfc);
 hold on
   prismatic_joint_axis(0.0015,[b_s6_o(1),b_s6_o(2),b_s6_o(3)],[p_s6_o(1),p_s6_o(2),p_s6_o(3)],jointfc,jointec,linkfc);
 hold on
 delta_1=p_s1_o-b_s1_o/norm(p_s1_o-b_s1_o);
 delta_2=p_s3_o-b_s3_o/norm(p_s3_o-b_s3_o);
 delta_3=p_s5_o-b_s5_o/norm(p_s5_o-b_s5_o);
 sphere_joint_axis(10,delta_1,0.004,eye(3),b_s1_o',jointfc,jointec);
sphere_joint_axis(10,delta_2,0.004,eye(3),b_s3_o',jointfc,jointec);
sphere_joint_axis(10,delta_3,0.004,eye(3),b_s5_o',jointfc,jointec);

 delta_11=p_s1_o-b_s1_o/norm(p_s1_o-b_s1_o);
 delta_22=p_s3_o-b_s3_o/norm(p_s3_o-b_s3_o);
 delta_33=p_s5_o-b_s5_o/norm(p_s5_o-b_s5_o);
 sphere_joint_axis(10,delta_11,0.004,eye(3),p_s1_o',jointfc,jointec);
sphere_joint_axis(10,delta_22,0.004,eye(3),p_s3_o',jointfc,jointec);
sphere_joint_axis(10,delta_33,0.004,eye(3),p_s5_o',jointfc,jointec);
hold on
%% plot plane
% % axes('nextplot', 'add');
% drawplanethreepoint(p_s1_o,p_s3_o,p_s5_o)
% hold on
%  plot3([b_s1_o(1) p_s1_o(1)]',[b_s1_o(2) p_s1_o(2)]',[b_s1_o(3) p_s1_o(3)]','Color' , 'k', 'LineWidth', 2.5, 'LineStyle','-');
%  hold on
%   plot3([b_s2_o(1) p_s2_o(1)]',[b_s2_o(2) p_s2_o(2)]',[b_s2_o(3) p_s2_o(3)]','Color' , 'k', 'LineWidth', 2.5, 'LineStyle','-');
%   hold on
%   plot3([b_s3_o(1) p_s3_o(1)]',[b_s3_o(2) p_s3_o(2)]',[b_s3_o(3) p_s3_o(3)]','Color' , 'k', 'LineWidth', 2.5, 'LineStyle','-');
%   hold on
%   plot3([b_s4_o(1) p_s4_o(1)]',[b_s4_o(2) p_s4_o(2)]',[b_s4_o(3) p_s4_o(3)]','Color' , 'k', 'LineWidth', 2.5, 'LineStyle','-');
%   hold on
%   plot3([b_s5_o(1) p_s5_o(1)]',[b_s5_o(2) p_s5_o(2)]',[b_s5_o(3) p_s5_o(3)]','Color' , 'k', 'LineWidth', 2.5, 'LineStyle','-');
%   hold on
%   plot3([b_s6_o(1) p_s6_o(1)]',[b_s6_o(2) p_s6_o(2)]',[b_s6_o(3) p_s6_o(3)]','Color' , 'k', 'LineWidth', 2.5, 'LineStyle','-');
%   hold on
  view(10.126508917753073,39.679461181786593)
%   drawnow
  axis off
  light('Position',[1 0 1],'Style','local')
  pause(0.1)
   hold off
%  hold on
%     currFrame= getframe(gcf);
%            writeVideo(vidObj, currFrame);
%       end
%      end
%  end
%      close(gcf)
%  close(vidObj);
p_os=0.001*[0;0;72];
b_s1_o=0.001*[0;129;0];b_s2_o=0.001*[0;129;0];b_s3_o=0.001*[-111.7;-64.5;0];b_s4_o=0.001*[-111.7;-64.5;0];
b_s5_o=0.001*[111.7;-64.5;0];b_s6_o=0.001*[111.7;-64.5;0];
p_s1_s=0.001*[42.4;24.5;51];p_s6_s=0.001*[42.4;24.5;51];p_s2_s=0.001*[-42.4;24.5;51];p_s3_s=0.001*[-42.4;24.5;51];
p_s4_s=0.001*[0;-49;51];p_s5_s=0.001*[0;-49;51];

%% jacobian or wrench matrix
% syms phi theta sigma real
% R_os=Rz(phi)*Ry(theta)*Rz(sigma-phi);
%% cable space jacobian
L_s1v= (mtimes(R_os,p_s1_s)-b_s1_o);
L_s2v= (mtimes(R_os,p_s2_s)-b_s2_o);
L_s3v= (mtimes(R_os,p_s3_s)-b_s3_o);
L_s4v= (mtimes(R_os,p_s4_s)-b_s4_o);
L_s5v= (mtimes(R_os,p_s5_s)-b_s5_o);
L_s6v= (mtimes(R_os,p_s6_s)-b_s6_o);

 u1=L_s1v/l_s1
 u2=L_s2v/l_s2;
 u3=L_s3v/l_s3;
 u4=L_s4v/l_s4;
 u5=L_s5v/l_s5;
 u6=L_s6v/l_s6;
 
 J_sho=[(cross(mtimes(R_os,p_s1_s),u1)).';
    (cross(mtimes(R_os,p_s2_s),u2)).';
    (cross(mtimes(R_os,p_s3_s),u3)).';
    (cross(mtimes(R_os,p_s4_s),u4)).';
    (cross(mtimes(R_os,p_s5_s),u5)).';
    (cross(mtimes(R_os,p_s6_s),u6)).';]
Wrench_sho=J_sho'
% A1=det(J_sho(1:3,1:3))
% A1=det(J_sho(4:6,1:3))
% 
% y=-(det(transpose(J_sho)*J_sho)); % manipulability

% %% n= 3 dof m=6 cable 3x6 matrix
% wld_position_base_pts=0.001*[b_s1_o,b_s2_o,b_s3_o,b_s4_o,b_s5_o,b_s6_o]'; % converting to meters from mm
% % tool_position_plat_pts=0.001*[p_s1_o, p_s3_o, p_s3_o, p_s4_o, p_s5_o, p_s6_o]'; % note here you have put transformed coordinates.
% tool_position_plat_pts=0.001*[p_s1_s, p_s3_s, p_s3_s, p_s4_s, p_s5_s, p_s6_s]'; % note here you have put  coordinates wrt S.
% 
% minimum_wrench = ([0.0; 0.0; 2.0 * 9.81]); %lets say platform must support it's weight of 2kg includes docking forces
% minimum_wrench_sub =[0.0; 0.0; 5.0 * 9.81];
%  % max and min cable limits
%     t_max = 120;  % maximum cable tension found form motors
%     t_min = 2 ; % need to ensure cable tension is above zero
%     tension_space_Vrep = tension_space_polytope(t_min,t_max,6);  % this defines a 'box' in tension space that has all feasible tensions
%    % An end effector position
%    
%       wld_pose_tool=eye(4);
% %     wld_pose_tool=Rz(pi/3)*Ry(pi/3)*Rz(0-pi/3);
%      wld_pose_tool(1:3,4)=0.001*[0,0,72];
%       W = wrench_matrixshoulder(wld_position_base_pts,wld_pose_tool, tool_position_plat_pts);  % Wt + we=0 
% %         figure
%  hold on
% %        polycart=polytope_Cartesian(Wrench_sho(1:3,:),tension_space_Vrep) %AWS available wrench set
% shift=center_plat';
%         polycart=polytope_Cartesian(Wrench_sho,tension_space_Vrep,0.01)+shift; 
% %      polycart=polytope_Cartesian(W(1:3,:),tension_space_Vrep) %AWS available wrench set
%          polycart.plot('wire',true)
%  light('Position',[1 0 1],'Style','local')
%  xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
%     polycart.plot('color','salmon','Alpha',0.3,'LineStyle','-','LineWidth',1.0,'Marker','o')
% 