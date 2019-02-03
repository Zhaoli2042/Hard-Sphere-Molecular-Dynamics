function a=UPPOS(n,pos,vel,tij,box,xx)
% update position after collision
for i=1:1:n
    pos(i,1)=pos(i,1)+vel(i,1)*tij;
    pos(i,2)=pos(i,2)+vel(i,2)*tij;
    % minimum image convention
    pos(i,1)=pos(i,1)-box*round(pos(i,1)/box);
    pos(i,2)=pos(i,2)-box*round(pos(i,2)/box);
end
a=pos;
% dump positions to txt file
ATOMS=linspace(1,n,n)';
id=ones(n,1);
for m = 1:1:n
    id(m)=id(m)-1+m;
end
matrix=[ATOMS id pos];
fid=fopen(sprintf('pos_%.0f.txt',xx),'w');
fprintf(fid,'ITEM: TIMESTEP\n');
fprintf(fid,'%.0f\n',1000);
fprintf(fid,'ITEM: NUMBER OF ATOMS\n');
fprintf(fid,'%.0f\n',n);
fprintf(fid,'ITEM: BOX BOUNDS pp pp pp\n');
% print box information to txt file
for x=1:1:3
    fprintf(fid,'%.0f %.0f\n',-0.5*box,0.5*box);
end
fprintf(fid,'ITEM: ATOMS id type xu yu zu \n');
% print position information to txt file
aa=size(matrix);aa=aa(1);
for x=1:1:aa
    fprintf(fid,'%.0f %.0f %.6f %.6f %.6f\n',matrix(x,1),matrix(x,2),matrix(x,3),matrix(x,4),0);
end