%the main vortex plotting routine, run this with a filenumber as an input, e.g.
%vortex_plot(1)
function vortex_boundary_plot(filenumber,varargin)
global dims box_size
global x y z
global f u u_mf v_curv
global ux uy uz
global number_of_particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('LineWidth', 1, @isscalar);
p.addParamValue('LineColor', 'k', @ischar);
p.addParamValue('LineStyle','-', @ischar);
p.addParamValue('MarkerColor','k', @ischar);
p.addParamValue('MaxPoints',30, @isscalar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxPoints=p.Results.MaxPoints;
rainbow=0;
if strcmp(p.Results.LineColor,'velocity')
  %scale velocity into a colormap
  store_caxis=([min(u(u>0)) max(u)]);
  u=u-min(u(u>0));
  rainbow_scale=199/max(u) ;
  u=u*rainbow_scale;
  rainbow_val=u;
  rainbowcmap=colormap(jet(200)); 
  rainbow=1;
elseif strcmp(p.Results.LineColor,'velocityx')
  %scale velocity into a colormap
  store_caxis=([min(ux) max(ux)]);
  ux=ux-min(ux);
  rainbow_scale=199/max(ux) ;
  ux=ux*rainbow_scale;
  rainbow_val=ux;
  rainbowcmap=colormap(jet(200)); 
  rainbow=1;
elseif strcmp(p.Results.LineColor,'velocityy')
  %scale velocity into a colormap
  store_caxis=([min(uy) max(uy)]);
  ux=ux-min(uy);
  rainbow_scale=199/max(uy) ;
  uy=uy*rainbow_scale;
  rainbow_val=uy;
  rainbowcmap=colormap(jet(200)); 
  rainbow=1;
elseif strcmp(p.Results.LineColor,'velocityz')
  %scale velocity into a colormap
  store_caxis=([min(uz) max(uz)]);
  uz=uz-min(uz);
  rainbow_scale=199/max(uz) ;
  uz=uz*rainbow_scale;
  rainbow_val=uz;
  rainbowcmap=colormap(jet(200)); 
  rainbow=1;
end
%deal with points on lower boundary
counter1=0;
for j=1:number_of_particles
  if abs(y(j)+dims(3)/2.)<0.01
   counter1=counter1+1 ;
   startingpos(counter1)=j;
  end
end 
%now plot loops
for i=1:counter1
  next=startingpos(i);
  plot=1;
  counter2=0;
  for j=1:number_of_particles
    nnext=f(next);
    dist=sqrt((x(next)-x(nnext))^2+(y(next)-y(nnext))^2+(z(next)-z(nnext))^2);
    if dist>3*dims(1)
     plot=0;
    end
    counter2=counter2+1;
    dummy_array(counter2,1)=x(next);
    dummy_array(counter2,2)=y(next);
    dummy_array(counter2,3)=z(next);
    if rainbow
      dummy_array(counter2,4)=rainbow_val(next);
    end
    if nnext==next
      break
    else
      next=nnext;
    end
  end
  if plot==1  && length(dummy_array)<MaxPoints && size(dummy_array,1)>3
  if rainbow
    for j=1:size(dummy_array,1)-1
    dummy_x(1,1)=dummy_array(j,1);
    dummy_x(2,1)=dummy_array(j+1,1);
    dummy_x(1,2)=dummy_array(j,2);
    dummy_x(2,2)=dummy_array(j+1,2);
    dummy_x(1,3)=dummy_array(j,3);
    dummy_x(2,3)=dummy_array(j+1,3);
    plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
    'Color',rainbowcmap(max(1,ceil(dummy_array(j,4))),:),'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)
    end
  else
    plot3(dummy_array(:,1),dummy_array(:,2),dummy_array(:,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
    'Color',p.Results.LineColor,'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)
  end
  end
  hold on
  clear dummy_array
end

%deal with points on upper boundary
counter1=0;
for j=1:number_of_particles
  if abs(y(j)-dims(3)/2.)<0.01
   counter1=counter1+1 ;
   startingpos(counter1)=j;
  end
end 
%now plot loops
for i=1:counter1
  next=startingpos(i);
  plot=1;
  counter2=0;
  for j=1:number_of_particles
    nnext=f(next);
    dist=sqrt((x(next)-x(nnext))^2+(y(next)-y(nnext))^2+(z(next)-z(nnext))^2);
    if dist>3*dims(1)
     plot=0;
    end
    counter2=counter2+1;
    dummy_array(counter2,1)=x(next);
    dummy_array(counter2,2)=y(next);
    dummy_array(counter2,3)=z(next);
    if rainbow
      dummy_array(counter2,4)=rainbow_val(next);
    end
    if nnext==next
      break
    else
      next=nnext;
    end
  end
  if plot==1  && length(dummy_array)<MaxPoints && size(dummy_array,1)>3
  if rainbow
    for j=1:size(dummy_array,1)-1
    dummy_x(1,1)=dummy_array(j,1);
    dummy_x(2,1)=dummy_array(j+1,1);
    dummy_x(1,2)=dummy_array(j,2);
    dummy_x(2,2)=dummy_array(j+1,2);
    dummy_x(1,3)=dummy_array(j,3);
    dummy_x(2,3)=dummy_array(j+1,3);
    plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
    'Color',rainbowcmap(max(1,ceil(dummy_array(j,4))),:),'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)
    end
  else
    plot3(dummy_array(:,1),dummy_array(:,2),dummy_array(:,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
    'Color',p.Results.LineColor,'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)
  end
  end
  hold on
  clear dummy_array
end
axis([-box_size(1)/2 box_size(1)/2 -box_size(2)/2 box_size(2)/2 -box_size(3)/ 2 box_size(3)/2])
  box on
axis square
  if rainbow
    caxis(store_caxis)
    colorbar
  end
camproj('perspective')
rotate3d on
daspect([1 1 1])
hold off

