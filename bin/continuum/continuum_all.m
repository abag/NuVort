function continuum_all(files,N)
if nargin==1
  N=16
end
global dims box_size
global Density Curv Friction allFriction Velocity allVelocity
global GDensity GCurv GFriction GallFriction GVelocity GallVelocity
global N
GDensity(1:N,1:N,1:N)=0.;
GCurv(1:N,1:N,1:N)=0.;
GFriction(1:N,1:N,1:N)=0.;
GallFriction(1:N,1:N,1:N,1:3)=0.;
GVelocity(1:N,1:N,1:N)=0.;
GallVelocity(1:N,1:N,1:N,1:3)=0.;
for i=files
  vortex_continuum(i,'N',N)
  GDensity=GDensity+Density;
  GCurv=GCurv+Curv ; 
  GFriction=GFriction+Friction;
  GallFriction=GallFriction+allFriction;
  GVelocity=GVelocity+Velocity;
  GallVelocity=GallVelocity+allVelocity;
end
GDensity=GDensity/length(files);
GCurv=GCurv/length(files);
GFriction=GFriction/length(files);
GallFriction=GallFriction/length(files);
GVelocity=GVelocity/length(files);
GallVelocity=GallVelocity/length(files);

return
%-------------plots---------------------
figure('Name','vortex density in xz-plane')
pcolor(xbox(1:N),zbox(1:N),squeeze(sum(GDensity,2))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','vortex density in xy-plane')
pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Density(:,:,:),3))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','vortex density in yz-plane')
pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Density(:,:,:),1))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','curvature in xz-plane')
pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Curv,2))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','curvature in xy-plane')
pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Curv(:,:,:),3))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','curvature in yz-plane')
pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Curv(:,:,:),1))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','mutual friction in xz-plane')
pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Friction,2))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','mutual friction in xy-plane')
pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Friction(:,:,:),3))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','mutual friction in yz-plane')
pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Friction(:,:,:),1))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);

figure('Name','mutual friction (x,y,z comp.) in xz-plane')
subplot(3,1,1)
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,1),2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
subplot(3,1,2)
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,2),2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
subplot(3,1,3)
  pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,3),2))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
   
figure('Name','mutual friction (x,y,z comp.) in xy-plane')
subplot(3,1,1)
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allFriction(:,:,:,1),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
subplot(3,1,2)
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allFriction(:,:,:,2),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
subplot(3,1,3)
  pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allFriction(:,:,:,3),3))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16);
   
figure('Name','mutual friction (x,y,z comp.) in yz-plane')
subplot(1,3,1)
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,1),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
subplot(1,3,2)
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,2),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
subplot(1,3,3)
  pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allFriction(:,:,:,3),1))') ; shading interp
  daspect([1 1 1]) ; colorbar
  xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16);
  
figure('Name','total velocity in xz-plane')
pcolor(xbox(1:N),zbox(1:N),squeeze(sum(Velocity,2))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)


figure('Name','total velocity in xy-plane')
pcolor(xbox(1:N),ybox(1:N),squeeze(sum(Velocity(:,:,:),3))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
  
figure('Name','total velocity in yz-plane')
pcolor(ybox(1:N),zbox(1:N),squeeze(sum(Velocity(:,:,:),1))') ; shading interp
daspect([1 1 1]) ; colorbar
xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)

figure('Name','velocity (x,y,z comp.) in xz-plane')
subplot(3,1,1)
 pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,1),2))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
subplot(3,1,2)
 pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,2),2))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
subplot(3,1,3)
 pcolor(xbox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,3),2))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('x','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)

figure('Name','velocity (x,y,z comp.) in xy-plane')
subplot(3,1,1)
 pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allVelocity(:,:,:,1),3))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
subplot(3,1,2)
 pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allVelocity(:,:,:,2),3))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)
subplot(3,1,3)
 pcolor(xbox(1:N),ybox(1:N),squeeze(sum(allVelocity(:,:,:,3),3))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('x','FontSize',16) ; ylabel('y','FontSize',16) ; set(gca,'FontSize',16)

figure('Name','velocity (x,y,z comp.) in yz-plane')
subplot(1,3,1)
 pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,1),1))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
subplot(1,3,2)
 pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,2),1))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)
subplot(1,3,3)
 pcolor(ybox(1:N),zbox(1:N),squeeze(sum(allVelocity(:,:,:,3),1))') ; shading interp
 daspect([1 1 1]) ; colorbar
 xlabel('y','FontSize',16) ; ylabel('z','FontSize',16) ; set(gca,'FontSize',16)