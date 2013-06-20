function mean_density(files)
global Density box_size
mean_d_array=0.;
for i=files
    N=vortex_density(i);
    mean_d_array= mean_d_array+Density;
end
mean_d_array=mean_d_array/length(files);
xbox=-box_size(1)/2:box_size(1)/N:box_size(1)/2;
ybox=-box_size(2)/2:box_size(2)/N:box_size(2)/2;
zbox=-box_size(3)/2:box_size(3)/N:box_size(3)/2;
pcolor(xbox(1:N),zbox(1:N),squeeze(sum(mean_d_array,2))) ; shading interp
daspect([1 1 1])
figure
pcolor(xbox(1:N),ybox(1:N),squeeze(sum(mean_d_array(:,:,:),3))) ; shading interp
daspect([1 1 1])
figure
pcolor(ybox(1:N),zbox(1:N),squeeze(sum(mean_d_array(:,:,:),1))) ; shading interp
daspect([1 1 1])