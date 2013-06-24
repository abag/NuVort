function N=vortex_density(filenumber,varargin)
global dims box_size
global x y z
global f u u2 v_curv u_mf
global u_mf_x u_mf_y u_mf_z
global ux uy uz
global number_of_particles
global Density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('N',16, @isscalar);
p.addParamValue('Plots','off', @ischar);
p.addParamValue('PlotPlane','all', @ischar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create a NxNxN array
N=p.Results.N;
xbox=-box_size(1)/2:box_size(1)/N:box_size(1)/2;
ybox=-box_size(2)/2:box_size(2)/N:box_size(2)/2;
zbox=-box_size(3)/2:box_size(3)/N:box_size(3)/2;
Density(1:N,1:N,1:N)=0.;
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
if strcmp(p.Results.Plots,'density')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
    figure('Name','vortex density in xz-plane')
    pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Density,2))') ; shading interp
    daspect([1 1 1]) ; colorbar
    xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
    figure('Name','vortex density in xy-plane')
    pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Density(:,:,:),3))') ; shading interp
    daspect([1 1 1]) ; colorbar
    xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
    figure('Name','vortex density in yz-plane')
    pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Density(:,:,:),1))') ; shading interp
    daspect([1 1 1]) ; colorbar
    xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
elseif  strcmp(p.Results.Plots,'curv')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure('Name','curvature in xz-plane')
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Curv,2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure('Name','curvature in xy-plane')
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Curv(:,:,:),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure('Name','curvature in yz-plane')
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Curv(:,:,:),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
elseif  strcmp(p.Results.Plots,'friction')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure('Name','mutual friction in xz-plane')
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Friction,2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure('Name','mutual friction in xy-plane')
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Friction(:,:,:),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure('Name','mutual friction in yz-plane')
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Friction(:,:,:),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
elseif  strcmp(p.Results.Plots,'allfriction')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure('Name','mutual friction (x,y,z comp.) in xz-plane')
  subplot(3,1,1)
   pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,1),2))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  subplot(3,1,2)
   pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,2),2))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  subplot(3,1,3)
   pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,3),2))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure('Name','mutual friction (x,y,z comp.) in xy-plane')
  subplot(3,1,1)
   pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allFriction(:,:,:,1),3))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
  subplot(3,1,2)
   pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allFriction(:,:,:,2),3))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
  subplot(3,1,3)
   pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allFriction(:,:,:,3),3))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure('Name','mutual friction (x,y,z comp.) in yz-plane')
  subplot(1,3,1)
   pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,1),1))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  subplot(1,3,2)
   pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,2),1))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  subplot(1,3,3)
   pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,3),1))') ; shading interp
   daspect([1 1 1]) ; colorbar
   xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  end
elseif  strcmp(p.Results.Plots,'velocity')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure('Name','total velocity in xz-plane')
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Velocity,2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure('Name','total velocity in xy-plane')
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Velocity(:,:,:),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure('Name','total velocity in yz-plane')
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Velocity(:,:,:),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  end
elseif  strcmp(p.Results.Plots,'allvelocity')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure('Name','velocity (x,y,z comp.) in xz-plane')
  subplot(3,1,1)
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,1),2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  subplot(3,1,2)
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,2),2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  subplot(3,1,3)
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,3),2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure('Name','velocity (x,y,z comp.) in xy-plane')
  subplot(3,1,1)
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allVelocity(:,:,:,1),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
  subplot(3,1,2)
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allVelocity(:,:,:,2),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
  subplot(3,1,3)
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allVelocity(:,:,:,3),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure('Name','velocity (x,y,z comp.) in yz-plane')
  subplot(1,3,1)
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,1),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  subplot(1,3,2)
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,2),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  subplot(1,3,3)
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,3),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
  end
end
