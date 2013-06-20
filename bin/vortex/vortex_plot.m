%the main vortex plotting routine, run this with a filenumber as an input, e.g.
%vortex_plot(1)
function vortex_plot(filenumber,varargin)
global dims box_size
global x y z
global f u u_mf v_curv
global number_of_particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('LineWidth', 1, @isscalar);
p.addParamValue('Views', 1, @isscalar);
p.addParamValue('LineColor', 'k', @ischar);
p.addParamValue('LineStyle','-', @ischar);
p.addParamValue('MarkerColor','k', @ischar);
p.addParamValue('OverHead','off', @ischar);
p.addParamValue('Annotate','off', @ischar);
parse(p,varargin{:});
if p.Results.Views~=[1 2 4]
  disp('Views must be set to either 1,2 or 4')
  return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vortex_load(filenumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
elseif strcmp(p.Results.LineColor,'friction')
  store_caxis=([min(u_mf(u_mf>0)) max(u_mf)]);
  u_mf=u_mf-min(u_mf(u_mf>0));
  rainbow_scale=199/max(u_mf) ;
  u_mf=u_mf*rainbow_scale;
  rainbow_val=u_mf;
  rainbowcmap=colormap(jet(200)); 
  rainbow=1;
elseif strcmp(p.Results.LineColor,'curvature')
  store_caxis=([min(v_curv(v_curv>0)) max(v_curv)]);
  v_curv=v_curv-min(v_curv(v_curv>0));
  rainbow_scale=199/max(v_curv) ;
  v_curv=v_curv*rainbow_scale;
  rainbow_val=v_curv;
  rainbowcmap=colormap(jet(200)); 
  rainbow=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for v=1:(p.Results.Views)
  subplot(floor(sqrt((p.Results.Views))),ceil(sqrt((p.Results.Views))),v)
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
        if rainbow
          plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
            'Color',rainbowcmap(max(1,ceil(rainbow_val(j))),:),'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)
        else
          plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),p.Results.LineStyle,'LineWidth',p.Results.LineWidth,...
            'Color',p.Results.LineColor,'MarkerFaceColor',p.Results.MarkerColor,'MarkerEdgeColor',p.Results.MarkerColor)
        end
        hold on
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %set axis
  axis([-box_size(1)/2 box_size(1)/2 -box_size(2)/2 box_size(2)/2 -box_size(3)/2 box_size(3)/2])
  box on
  if rainbow
    caxis(store_caxis)
    colorbar
  end
  if (v==1) ; view(-45,30) ; end
  if (v==2) ; view(0,90) ; end
  if (v==3) ; view(90,0) ; end
  if (v==4) ; view(45,30) ; end
  if strcmp(p.Results.OverHead,'on')
    view(0,90)
  end
  set(gca,'FontSize',16)
  if strcmp(p.Results.Annotate,'off')
    set(gca,'xtick',[]) ; set(gca,'ytick',[]) ; set(gca,'ztick',[])
  end
  axis square
  camproj('perspective')
  rotate3d on
  daspect([1 1 1])
  hold off
end
hold off
clear all
