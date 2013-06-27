function continuum_plot2D(varargin)
global dims box_size
global GDensity GCurv GFriction GallFriction GVelocity GallVelocity
global N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('Plots','off', @ischar);
p.addParamValue('PlotPlane','xy', @ischar);
p.addParamValue('PlotStyle','mean', @ischar);
p.addParamValue('SnapPos',N/2, @isscalar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbox=-box_size(1)/2:box_size(1)/N:box_size(1)/2;
ybox=-box_size(2)/2:box_size(2)/N:box_size(2)/2;
zbox=-box_size(3)/2:box_size(3)/N:box_size(3)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First what do we want to plot?
if strcmp(p.Results.Plots,'density')
  plot_array=GDensity ;
elseif strcmp(p.Results.Plots,'curv')
  plot_array=GCurv ;
elseif strcmp(p.Results.Plots,'friction')
  plot_array=GFriction ;
elseif strcmp(p.Results.Plots,'frictionx')
  plot_array=GallFriction(:,:,:,1) ;
elseif strcmp(p.Results.Plots,'frictiony')
  plot_array=GallFriction(:,:,:,2) ;
elseif strcmp(p.Results.Plots,'frictionz')
  plot_array=GallFriction(:,:,:,3) ;
elseif strcmp(p.Results.Plots,'velocity')
  plot_array=GVelocity ;
elseif strcmp(p.Results.Plots,'velocityx')
  plot_array=GallVelocity(:,:,:,1) ;
elseif strcmp(p.Results.Plots,'velocityy')
  plot_array=GallVelocity(:,:,:,2) ;
elseif strcmp(p.Results.Plots,'velocityz')
  plot_array=GallVelocity(:,:,:,3) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now what is the plane we want to plot in?
SnapPos=p.Results.SnapPos ;
if strcmp(p.Results.PlotPlane,'xy')
  if strcmp(p.Results.PlotStyle,'mean')
    plot_array2=squeeze(mean(plot_array(:,:,:),3))';
  elseif strcmp(p.Results.PlotStyle,'snap')
    plot_array2=squeeze(plot_array(:,:,SnapPos))' ;
  end
  x1=xbox ; x2=ybox; x3=zbox ;
  xlab='x' ; ylab='y' ; zlab='z' ;
elseif strcmp(p.Results.PlotPlane,'xz')
  if strcmp(p.Results.PlotStyle,'mean')
    plot_array2=squeeze(mean(plot_array(:,:,:),2))';
  elseif strcmp(p.Results.PlotStyle,'snap')
    plot_array2=squeeze(plot_array(:,SnapPos,:))' ;
  end
  x1=xbox ; x2=zbox; x3=ybox ;
  xlab='x' ; ylab='z' ; zlab='y' ;
elseif strcmp(p.Results.PlotPlane,'yz')
  if strcmp(p.Results.PlotStyle,'mean')
    plot_array2=squeeze(mean(plot_array(:,:,:),1))';
  elseif strcmp(p.Results.PlotStyle,'snap')
    plot_array2=squeeze(plot_array(SnapPos,:,:))' ;
  end
  x1=ybox ; x2=zbox; x3=xbox ;
  xlab='y' ; ylab='z' ; zlab='x' ;
end
if strcmp(p.Results.PlotStyle,'mean')
  titlename=['mean ', p.Results.Plots,' plotted in ',p.Results.PlotPlane,' plane']
elseif strcmp(p.Results.PlotStyle,'snap');
  pos=num2str(x3(SnapPos+1)) ;
  titlename=[p.Results.Plots,' plotted in ',p.Results.PlotPlane,' plane at ',zlab,'=',pos];
end
figure('Name',titlename)
pcolor(x1(1:N),x2(1:N),plot_array2) ; shading interp
daspect([1 1 1])
colorbar
xlabel(xlab,'FontSize',16)
ylabel(ylab,'FontSize',16)
set(gca,'FontSize',16)




