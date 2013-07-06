function [L_inside L_outside L_total vort_rms]=kdtree_vort_smooth(filenumber,plot,eval_rms)
global dims
global x y z
global f
global number_of_particles
global vort total_length vort_rms
if nargin<2
  plot=1;
  eval_rms=2.;
elseif nargin<3
  eval_rms=2.;
end
%global test
vortex_load(filenumber);
counter=1;
for j=1:number_of_particles
  if f(j)==0 ; continue ; end
  dist=sqrt((x(j)-x(round(f(j))))^2+(y(j)-y(round(f(j))))^2+(z(j)-z(round(f(j))))^2);
  if (dist>0.5*min(dims(2:4))) ; continue ; end
  varray(counter,1)=x(j);
  varray(counter,2)=y(j);
  varray(counter,3)=z(j);
  varray(counter,4)=x(round(f(j)))-x(j);
  varray(counter,5)=y(round(f(j)))-y(j);
  varray(counter,6)=z(round(f(j)))-z(j);
  varray(counter,7:9)=0.;
  counter=counter+1;
end
l_varray=length(varray);
%----------------------------------------------------------------------
ns = createns(varray(:,1:3),'nsmethod','kdtree');
%now we need to define the mesh
%let's get the smoothing length
LL=sum(sqrt(varray(:,4).^2+varray(:,5).^2+varray(:,6).^2))
delta=sqrt((dims(2)*dims(3)*dims(4))/LL);
for i=1:l_varray ;
  [ind ddist]=knnsearch(ns,varray(i,1:3),'k',800);
  if max(ddist)<2*delta
    disp('i think more points are needed in KD tree')
  end
  for l=1:length(ind)
    factor=SPH_kernel_new(ddist(l),delta);
    varray(i,7)=varray(i,7)+factor*varray(ind(l),4);
    varray(i,8)=varray(i,8)+factor*varray(ind(l),5);
    varray(i,9)=varray(i,9)+factor*varray(ind(l),6);
  end
  clear ind ddist
end ; 
vort=sqrt(varray(:,7).^2+varray(:,8).^2+varray(:,9).^2);
vort_rms=sqrt(sum(varray(:,7).^2+varray(:,8).^2+varray(:,9).^2)/l_varray);
total_length=sqrt(varray(:,4).^2+varray(:,5).^2+varray(:,6).^2);
vort_rms
  
  vort_norm=vort-min(vort);
  rainbow_scale=299/max(vort_norm);
  vort_norm=vort_norm*rainbow_scale;
  store_caxis=([min(vort) max(vort)]);
  rainbowcmap=colormap(hot(400));
  for i=1:l_varray
    plot3([varray(i,1) (varray(i,4)+varray(i,1))],[varray(i,2) (varray(i,5)+varray(i,2))],[varray(i,3) (varray(i,6)+varray(i,3))],'-','Color',rainbowcmap(max(1,ceil(vort_norm(i))),:),'LineWidth',1.)
%     subplot(2,2,2)
%     if vort(i)>1*vort_rms
%           plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
%           [varray(i,2) (varray(i,5)+varray(i,2))],...
%           [varray(i,3) (varray(i,6)+varray(i,3))],...
%           '-','Color',rainbowcmap(max(1,ceil(10*vort(i))),:),...
%           'LineWidth',1.)
% 
%     hold on
%     end
%     subplot(2,2,4)
%     if vort(i)<1*vort_rms
%           plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
%           [varray(i,2) (varray(i,5)+varray(i,2))],...
%           [varray(i,3) (varray(i,6)+varray(i,3))],...
%           '-','Color',rainbowcmap(max(1,ceil(10*vort(i))),:),...
%           'LineWidth',1.)
% 
%     hold on
%     end
%     subplot(2,2,[1 3])
%     plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
%           [varray(i,2) (varray(i,5)+varray(i,2))],...
%           [varray(i,3) (varray(i,6)+varray(i,3))],...
%           '-','Color',rainbowcmap(max(1,ceil(10*vort(i))),:),...
%           'LineWidth',1.)

    hold on
  end
%   subplot(2,2,2)
%   camproj('perspective')
%   box on
%   axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
%   daspect([1 dims(7) dims(7)])
%   set(gca,'xtick',[])
%   set(gca,'ytick',[])
%   set(gca,'ztick',[])
%   subplot(2,2,4)
%   camproj('perspective')
%   box on
%   axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
%   daspect([1 dims(7) dims(7)])
%   set(gca,'xtick',[])
%   set(gca,'ytick',[])
%   set(gca,'ztick',[])
%   subplot(2,2,[1 3])
  box_size=dims(2:4)
  axis([-box_size(1)/2 box_size(1)/2 -box_size(2)/2 box_size(2)/2 -box_size(3)/ 2 box_size(3)/2])
  camproj('perspective')
  box on
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  caxis([0 1])
  colorbar
  axis square
  rotate3d on
  daspect([1 1 1])
  if 1==1
    whitebg('k')
    set(gcf,'InvertHardcopy','off');
  end
  hold off
end
function factor=SPH_kernel(r,h)
  q=r/(h);
  if q<=1 
    factor=(2-q)^3-4*(1-q)^3;
  elseif q<=2
    factor=(2-q)^3;
  else
    factor=0.;
  end
  factor=9.97E-4*factor*0.0163/(h^3);
end
function factor=SPH_kernel_new(r,h)
  %q=r/(2*h);
  q=r/h;
  if q<=1.
    factor=1-1.5*(q^2)+0.75*(q^3);
  elseif q<=2
    factor=0.25*((2-q)^3);
  else
    factor=0.;
  end
  factor=9.97E-4*factor/((2*pi*h^2)^(3/2));
end
function factor=gaussian_kernel(r,h)
  factor=exp(-r^2/(2*h^2));
  factor=9.97E-4*factor/((2*pi*h^2)^(3/2));
end
function factor=SPH_kernel_springel(r,h)
  q=r/(2*h);
  if q<=0.5 
    factor=1-6*(q^2)+6*(q^3);
  elseif q<=1
    factor=2*((1-q)^3);
  else
    factor=0.;
  end
  factor=9.97E-4*2.54647909*factor/(h^3);
end
