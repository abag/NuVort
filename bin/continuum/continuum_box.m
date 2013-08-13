function Density=continuum_box(filenumber,varargin)
global dims box_size
global x y z
global f u u2 v_curv u_mf
global u_mf_x u_mf_y u_mf_z
global ux uy uz
global number_of_particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('D',1., @isscalar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create a NxNxN array
D=p.Results.D;
Density=0.;
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
        if x(j)>-D/2. && x(j)<D/2. && ...
           y(j)>-D/2. && y(j)<D/2. && ...
           z(j)>-D/2. && z(j)<D/2.
           Density=Density+dist;
        end
      end
    end
end
%------------------------
Density=Density/(dims(2)*dims(3)*dims(4));
%------------------------

