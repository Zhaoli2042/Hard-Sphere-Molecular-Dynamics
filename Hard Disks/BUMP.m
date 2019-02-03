function [a,b]=BUMP(sigma,i,j,pos,vel,box)
% this shit updates velocities and virial
sigsq=sigma*sigma;
rxij=pos(i,1)-pos(j,1);
ryij=pos(i,2)-pos(j,2);
% minimum image convention
rxij=rxij-box*round(rxij/box);
ryij=ryij-box*round(ryij/box);
vxij=vel(i,1)-vel(j,1);
vyij=vel(i,2)-vel(j,2);
bij=rxij*vxij+ryij*vyij;
factor=bij/sigsq;
% Allen Tildesley eqn3.69
delvx=-factor*rxij;
delvy=-factor*ryij;
% A and T: eqn3.68
vel(i,1)=vel(i,1)+delvx;
vel(i,2)=vel(i,2)+delvy;
vel(j,1)=vel(j,1)-delvx;
vel(j,2)=vel(j,2)-delvy;
% update virial
w=delvx*rxij+delvy*ryij;
a=vel;
b=w;