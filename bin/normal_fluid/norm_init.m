function norm_init
close all
filename='data/norm_init_mesh.dat';
load data/nfm_dims.log;
msize=nfm_dims(1)
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
x=fread(fid,msize,'float64');
ux=fread(fid,msize^3,'float64');
disp('mean of ux')
mean(ux)
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
fclose(fid)
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
unorm=max(sqrt(ux.^2+uy.^2+uz.^2));
unorm=max(max(max(sqrt(ux.^2+uy.^2+uz.^2))));
figure('Name','u(x)')
  plot(squeeze(ux(:,32,32)))
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize,'initial normal fluid')
%plot slices of divergence of field
ux=permute(ux,[2 3 1]);
uy=permute(uy,[2 3 1]);
uz=permute(uz,[2 3 1]);
figure('Name','div(u)')
div = divergence(ux,uy,uz);
pcolor(squeeze(div(32,:,:))) ; shading interp
colorbar
%plot slices of curl of field
figure('Name','curl(u)')
curlz = curl(ux,uy,uz);
pcolor(squeeze(curlz(32,:,:))) ; shading interp
colorbar
[shearx1 shearx2 shearx3]=gradient(ux,x(2)-x(1));
[sheary1 sheary2 sheary3]=gradient(uy,x(2)-x(1));
[shearz1 shearz2 shearz3]=gradient(uz,x(2)-x(1));
Shear=abs(shearx3+shearz1);
figure
  p=patch(isosurface(x,x,x,Shear,0.999*max(max(max(Shear)))));
  isonormals(x,x,x,Shear, p)
  set(p, 'FaceColor', rgb('Green'), 'EdgeColor', 'none');
  alpha(0.4)
  daspect([1 1 1]); axis tight;
  camup([0 0 1 ]); campos([0.7686    0.1432    0.3043])
  camlight; lighting phong
  