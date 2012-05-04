%the main vortex plotting routine, run this with a filenumber as an input, e.g.
%vortex_plot(1)
function vortex_plot(filenumber,varargin)
global dims
global x y z
global f u u2
global number_of_particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('LineWidth', 1, @isscalar); 
p.addParamValue('LineColor', 'k', @ischar);
p.addParamValue('LineStyle','-', @ischar);
p.addParamValue('MarkerColor','k', @ischar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if (dist<0.5*dims(2))
      plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
            'Color',p.Results.LineColor,'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)  
    end
    hold on
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set axis
axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]) 
box on
set(gca,'FontSize',16)
axis square
camproj('perspective')
rotate3d on
hold off