function vortex_continuum(filenumber,varargin)
global dims box_size
global x y z
global f u u2 v_curv u_mf
global u_mf_x u_mf_y u_mf_z
global ux uy uz
global number_of_particles
global Density allLength Curv Friction allFriction Velocity allVelocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('N',16, @isscalar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create a NxNxN array
N=p.Results.N;
xbox=-box_size(1)/2:box_size(1)/N:box_size(1)/2;
ybox=-box_size(2)/2:box_size(2)/N:box_size(2)/2;
zbox=-box_size(3)/2:box_size(3)/N:box_size(3)/2;
Density(1:N,1:N,1:N)=0.;
allLength(1:N,1:N,1:N,1:3)=0.;
Curv(1:N,1:N,1:N)=0.;
Count(1:N,1:N,1:N)=0.;
Friction(1:N,1:N,1:N)=0.;
allFriction(1:N,1:N,1:N,1:3)=0.;
Velocity(1:N,1:N,1:N)=0.;
allVelocity(1:N,1:N,1:N,1:3)=0.;
for j=1:number_of_particles
    if round(f(j))==0
    else
      dummy_x(1,1)=x(j);
      dummy_x(2,1)=x(round(f(j)));
      dummy_x(1,2)=y(j);
      dummy_x(2,2)=y(round(f(j)));
      dummy_x(1,3)=z(j);
      dummy_x(2,3)=z(round(f(j)));
      dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
      if (dist<0.5*min(box_size))
        for ii=1:N
            for jj=1:N
                for kk=1:N
                    if x(j)>xbox(ii) && x(j)<xbox(ii+1) && ...
                       y(j)>ybox(jj) && y(j)<ybox(jj+1) && ...
                       z(j)>zbox(kk) && z(j)<zbox(kk+1)
                       Density(ii,jj,kk)=Density(ii,jj,kk)+dist;
                       allLength(ii,jj,kk,1)=allLength(ii,jj,kk,1)+abs(dummy_x(2,1)-dummy_x(1,1));
                       allLength(ii,jj,kk,2)=allLength(ii,jj,kk,2)+abs(dummy_x(2,2)-dummy_x(1,2));
                       allLength(ii,jj,kk,3)=allLength(ii,jj,kk,3)+abs(dummy_x(2,3)-dummy_x(1,3));
                       Count(ii,jj,kk)=Count(ii,jj,kk)+1;
                       Curv(ii,jj,kk)=Curv(ii,jj,kk)+v_curv(j);
                       Friction(ii,jj,kk)=Friction(ii,jj,kk)+u_mf(j);
                       allFriction(ii,jj,kk,1)=allFriction(ii,jj,kk,1)+u_mf_x(j);
                       allFriction(ii,jj,kk,2)=allFriction(ii,jj,kk,2)+u_mf_y(j);
                       allFriction(ii,jj,kk,3)=allFriction(ii,jj,kk,3)+u_mf_z(j);
                       Velocity(ii,jj,kk)=Velocity(ii,jj,kk)+u(j);
                       allVelocity(ii,jj,kk,1)=allVelocity(ii,jj,kk,1)+ux(j);
                       allVelocity(ii,jj,kk,2)=allVelocity(ii,jj,kk,2)+uy(j);
                       allVelocity(ii,jj,kk,3)=allVelocity(ii,jj,kk,3)+uz(j);
                    end
                end
            end
        end
      end
    end
end
%------------------------
allLength(:,:,:,1)=allLength(:,:,:,1)./Density;
allLength(:,:,:,2)=allLength(:,:,:,2)./Density;
allLength(:,:,:,3)=allLength(:,:,:,3)./Density;
tf=isnan(allLength);
allLength(tf)=0.;
%------------------------
Density=Density/(dims(2)*dims(3)*dims(4));
%------------------------
Curv=Curv./Count;
tf=isnan(Curv);
Curv(tf)=0.;
%------------------------
Friction=Friction./Count;
tf=isnan(Friction);
Friction(tf)=0.;
%------------------------
allFriction(:,:,:,1)=allFriction(:,:,:,1)./Count;
allFriction(:,:,:,2)=allFriction(:,:,:,2)./Count;
allFriction(:,:,:,3)=allFriction(:,:,:,3)./Count;
tf=isnan(allFriction);
allFriction(tf)=0.;
%------------------------
Velocity=Velocity./Count;
tf=isnan(Velocity);
Velocity(tf)=0.;
%------------------------
allVelocity(:,:,:,1)=allVelocity(:,:,:,1)./Count;
allVelocity(:,:,:,2)=allVelocity(:,:,:,2)./Count;
allVelocity(:,:,:,3)=allVelocity(:,:,:,3)./Count;
tf=isnan(allVelocity);
allVelocity(tf)=0.;
%------------------------

