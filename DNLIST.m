% Find collisions--update collisions
function [a,b]=DNLIST(sigma,ii,jj,pos,vel,coltim,partner,box,n,timbig)
sigsq=sigma^2;
for i=1:1:n
    if ((i == ii) || (i == jj) || (partner(i) == ii) ||(partner(i) == jj))
        coltim(i)=timbig; % set the collided times back to timebig
        for j=1:1:n
            if (j ~= i)
                rxj=pos(j,1);
                ryj=pos(j,2);
                rzj=pos(j,3);
                vxj=vel(j,1);
                vyj=vel(j,2);
                vzj=vel(j,3);
                rxij=pos(i,1)-rxj;
                ryij=pos(i,2)-ryj;
                rzij=pos(i,3)-rzj;
                % minimum image convention
                rxij=rxij-box*round(rxij/box);
                ryij=ryij-box*round(ryij/box);
                rzij=rzij-box*round(rzij/box);
                vxij=vel(i,1)-vxj;
                vyij=vel(i,2)-vyj;
                vzij=vel(i,3)-vzj;
                bij=rxij*vxij+ryij*vyij+rzij*vzij;
                if (bij < 0)
                    % only when bij <0 can yield positive time values
                    rijsq=rxij*rxij+ryij*ryij+rzij*rzij;
                    vijsq=vxij*vxij+vyij*vyij+vzij*vzij;
                    discr=bij*bij-vijsq*(rijsq-sigsq);
                    % at this case, only discr > 0 yields 2 positive roots
                    if (discr > 0)
                        tij=(-bij-sqrt(discr))/vijsq;
                        % use the calculated tij values to update collision times
                        if (tij < coltim(i))
                            coltim(i)=tij;
                            partner(i)=j; % set j as collision partner
                        end
                        if (tij < coltim(j))
                            coltim(j)=tij;
                            partner(j)=i;
                        end
                    end
                end
            end
        end
    end
end
a=coltim;
b=partner;
