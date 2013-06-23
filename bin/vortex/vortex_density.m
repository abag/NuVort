function N=vortex_density(filenumber,varargin)
global dims box_size
global x y z
global f u u2 v_curv u_mf
global number_of_particles
global Density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('Plots','off', @ischar);
p.addParamValue('PlotPlane','all', @ischar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create a NxNxN array
N=16;
xbox=-box_size(1)/2:box_size(1)/N:box_size(1)/2;
ybox=-box_size(2)/2:box_size(2)/N:box_size(2)/2;
zbox=-box_size(3)/2:box_size(3)/N:box_size(3)/2;
Density(N,N,N)=0.;
Curv(N,N,N)=0.;
Count(N,N,N)=0.;
Friction(N,N,N)=0.;
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
                    end
                end
            end
        end
      end
    end
end
Curv=Curv./Count;
tf=isnan(Curv);
Curv(tf)=0.;
Friction=Friction./Count;
tf=isnan(Friction);
Friction(tf)=0.;
if strcmp(p.Results.Plots,'density')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
    pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Density,2))) ; shading interp
    daspect([1 1 1])
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
    figure
    pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Density(:,:,:),3))) ; shading interp
    daspect([1 1 1])
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
    figure
    pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Density(:,:,:),1))) ; shading interp
    daspect([1 1 1])
  end
elseif  strcmp(p.Results.Plots,'curv')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Curv,2))) ; shading interp
  daspect([1 1 1])
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Curv(:,:,:),3))) ; shading interp
  daspect([1 1 1])
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Curv(:,:,:),1))) ; shading interp
  daspect([1 1 1])
  end
elseif  strcmp(p.Results.Plots,'friction')
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xz')
  figure
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Friction,2))) ; shading interp
  daspect([1 1 1])
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'xy')
  figure
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Friction(:,:,:),3))) ; shading interp
  daspect([1 1 1])
  end
  if strcmp(p.Results.PlotPlane,'all') || strcmp(p.Results.PlotPlane,'yz')
  figure
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Friction(:,:,:),1))) ; shading interp
  daspect([1 1 1])
  end
end
