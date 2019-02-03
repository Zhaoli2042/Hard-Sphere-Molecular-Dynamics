% Simple Hard Sphere Molecular Dynamics
clear all
clc
n=50;
rho=0.5;
vol=n/rho;
box=vol^(1/3);
sigma=1;  % This is the radius of particles
time=1000;
timestep=1;
% Generate random position and velocity
pos=zeros(n,3);
vel=zeros(n,3);
% place first particle on the shell
pos(1,:)=[0.5*box 0 0];
for i=2:1:n
    fprintf('Placing Atom Number %.0f\n',i);
    % minimum image convention
    isafe=0;
    while (isafe ~= i-1)
        isafe=0;
        rxi=0.5*box*(2*rand-1);
        ryi=0.5*box*(2*rand-1);
        rzi=0.5*box*(2*rand-1);
        rxi=rxi-box*round(rxi/box);
        ryi=ryi-box*round(ryi/box);
        rzi=rzi-box*round(rzi/box);
        for j=1:1:i-1
            % check overlap
            rxj=pos(j,1);
            ryj=pos(j,2);
            rzj=pos(j,3);
            rxij=rxi-rxj;
            ryij=ryi-ryj;
            rzij=rzi-rzj;
            % minimum image convention
            rxij=rxij-box*round(rxij/box);
            ryij=ryij-box*round(ryij/box);
            rzij=rzij-box*round(rzij/box);
            rijsq=rxij^2+ryij^2+rzij^2;
            if (rijsq > (sigma^2))
                isafe=isafe+1;
            end
        end
    end
    pos(i,:)=[rxi ryi rzi];
end
ATOMS=linspace(1,n,n)';
id=ones(n,1);
matrix=[ATOMS id pos];
fid=fopen('pos_test.txt','w');
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
    fprintf(fid,'%.0f %.0f %.6f %.6f %.6f\n',matrix(x,1),matrix(x,2),matrix(x,3),matrix(x,4),matrix(x,5));
end

% random velocity
for x=1:1:n
    for y=1:1:3
    vel(x,y)=2*rand-1;
    end
end
%% Here begins the Hard Sphere MD code
timbig=10^10; 
sigsq=sigma^2;
partner=zeros(n,1);
% initial calculations of collisions
rxi=pos(i,1); ryi=pos(i,2); rzi=pos(i,3);
vxi=vel(i,1); vyi=vel(i,2); vzi=vel(i,3);
for i=1:1:n
    coltim(i)=timbig;
end
for i=1:1:n-1
    for j=i+1:1:n
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
        if (bij < 0)
            rijsq=rxij*rxij+ryij*ryij+rzij*rzij;
            vijsq=vxij*vxij+vyij*vyij+vzij*vzij;
            discr=bij^2-vijsq*(rijsq-sigsq);
            if (discr > 0)
                tij = (-bij-sqrt(discr))/vijsq;
                if (tij < coltim(i))
                    coltim(i)=tij;
                    partner(i)=j;
                end
                if (tij < coltim(j))
                    coltim(j)=tij;
                    partner(j)=i;
                end
            end
        end
    end
end
% lets calculate 10 collisions
col=0;
tott=0;
for xx=1:1:50
% then locate the nearest collision
col=col+1;
tij=timbig;
for k=1:1:n
    if ((partner(k) ~= 0) && (coltim(k) < tij))
        tij=coltim(k);
        i=k;
    end
end
j=partner(i);
if (sum(partner) == 0)
    fprintf('No collision actually!');
    break
end
% now we have tij, we go update positions and velocities
coltim=coltim-tij; % after the first collision, time elapsed is tij
pos=UPPOS(n,pos,vel,tij,box,xx);
tott=tott+tij
% update new arrays

% then update velocity
[vel,w]=BUMP(sigma,i,j,pos,vel,box);
ii=i;
jj=j;
[coltim,partner]=DNLIST(sigma,ii,jj,pos,vel,coltim,partner,box,n,timbig);
end