function [a,b]=BUMP(sigma,i,j,pos,vel,box)
% this shit updates velocities and virial
sigsq=sigma*sigma;
rxij=pos(i,1)-pos(j,1);
ryij=pos(i,2)-pos(j,2);
rzij=pos(i,3)-pos(j,3);
% minimum image convention
rxij=rxij-box*round(rxij/box);
ryij=ryij-box*round(ryij/box);
rzij=rzij-box*round(rzij/box);
vxij=vel(i,1)-vel(j,1);
vyij=vel(i,2)-vel(j,2);
vzij=vel(i,3)-vel(j,3);
bij=rxij*vxij+ryij*vyij+rzij*vzij;
factor=bij/sigsq;
% Allen Tildesley eqn3.69
delvx=-factor*rxij;
delvy=-factor*ryij;
delvz=-factor*rzij;
% A and T: eqn3.68
vel(i,1)=vel(i,1)+delvx;
vel(i,2)=vel(i,2)+delvy;
vel(i,3)=vel(i,3)+delvz;
vel(j,1)=vel(j,1)-delvx;
vel(j,2)=vel(j,2)-delvy;
vel(j,3)=vel(j,3)-delvz;
% update virial
w=delvx*rxij+delvy*ryij+delvz*rzij;
a=vel;
b=w;