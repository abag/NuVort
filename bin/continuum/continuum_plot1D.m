function [x1 plot_array2]=continuum_plot1D(varargin)
global dims box_size
global GDensity GCurv GFriction GallFriction GVelocity GallVelocity GallLength
global N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('Plots','off', @ischar);
p.addParamValue('File','off', @ischar);
p.addParamValue('PlotPlane','x', @ischar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.Results.File,'off')
    disp('using array from memory')
else
    fin=[p.Results.File,'.mat'];
    load(fin)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbox=-box_size(1)/2:box_size(1)/N:box_size(1)/2;
ybox=-box_size(2)/2:box_size(2)/N:box_size(2)/2;
zbox=-box_size(3)/2:box_size(3)/N:box_size(3)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First what do we want to plot?
if strcmp(p.Results.Plots,'Lx')
  plot_array=GallLength(:,:,:,1) ;
  ylab='l_x/L' ;
elseif strcmp(p.Results.Plots,'Ly')
  plot_array=GallLength(:,:,:,2) ;
  ylab='l_y/L' ;
elseif strcmp(p.Results.Plots,'Lz')
  plot_array=GallLength(:,:,:,3) ;
  ylab='l_z/L' ;
elseif strcmp(p.Results.Plots,'density')
  plot_array=GDensity ;
  ylab='L' ;
elseif strcmp(p.Results.Plots,'curv')
  plot_array=GCurv ;
  ylab='\kappa' ;
elseif strcmp(p.Results.Plots,'friction')
  plot_array=GFriction ;
  ylab='|u_{mf}|' ;
elseif strcmp(p.Results.Plots,'frictionx')
  plot_array=GallFriction(:,:,:,1) ;
  ylab='u_{mf,x}' ;
elseif strcmp(p.Results.Plots,'frictiony')
  plot_array=GallFriction(:,:,:,2) ;
  ylab='u_{mf,y}' ;
elseif strcmp(p.Results.Plots,'frictionz')
  plot_array=GallFriction(:,:,:,3) ;
  ylab='u_{mf,z}' ;
elseif strcmp(p.Results.Plots,'velocity')
  plot_array=GVelocity ;
  ylab='|u|' ;
elseif strcmp(p.Results.Plots,'velocityx')
  plot_array=GallVelocity(:,:,:,1) ;
  ylab='u_x' ;
elseif strcmp(p.Results.Plots,'velocityy')
  plot_array=GallVelocity(:,:,:,2) ;
  ylab='u_y' ;
elseif strcmp(p.Results.Plots,'velocityz')
  plot_array=GallVelocity(:,:,:,3) ;
  ylab='u_z' ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now what is the plane we want to plot in?
if strcmp(p.Results.PlotPlane,'x')
  plot_array2=mean(squeeze(mean(plot_array(:,:,:),3)),2);
  %error_array2=std2(squeeze(std2(plot_array(:,:,:),3)),2);
  x1=xbox ; 
  xlab='x' ; 
elseif strcmp(p.Results.PlotPlane,'y')
  plot_array2=mean(squeeze(mean(plot_array(:,:,:),3)),1);
  %error_array2=std2(squeeze(std2(plot_array(:,:,:),3)),1);
  x1=ybox ;
  xlab='y' ; 
elseif strcmp(p.Results.PlotPlane,'z')
  plot_array2=mean(squeeze(mean(plot_array(:,:,:),1)),1);
  %error_array2=std2(squeeze(std2(plot_array(:,:,:),1)),1);
  x1=zbox ;
  xlab='z' ; 
end
titlename=['mean ', p.Results.Plots,' plotted along ',p.Results.PlotPlane,' axis']; 
figure('Name',titlename)
plot(x1(1:N),plot_array2,'LineWidth',2,'Color','k')
%errorbar(x1(1:N),plot_array2,error_array2,'LineWidth',2,'Color','k')
xlabel(xlab,'FontSize',16)
ylabel(ylab,'FontSize',16)
set(gca,'FontSize',16)




