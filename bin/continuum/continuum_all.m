function continuum_all(files,varargin)
global dims box_size
global Density Curv Friction allFriction Velocity allVelocity allLength
global GDensity GCurv GFriction GallFriction GVelocity GallVelocity GallLength
global N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validLogical = {'true','false'};
checkLogical = @(x) any(validatestring(x,validLogical));
p.addParamValue('N',16, @isscalar);
p.addParamValue('File','off', @ischar);
parse(p,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=p.Results.N;
GDensity(1:N,1:N,1:N)=0.;
GallLength(1:N,1:N,1:N,1:3)=0.;
GCurv(1:N,1:N,1:N)=0.;
GFriction(1:N,1:N,1:N)=0.;
GallFriction(1:N,1:N,1:N,1:3)=0.;
GVelocity(1:N,1:N,1:N)=0.;
GallVelocity(1:N,1:N,1:N,1:3)=0.;
for i=files
  vortex_continuum(i,'N',N)
  GDensity=GDensity+Density;
  GallLength=GallLength+allLength;
  GCurv=GCurv+Curv ; 
  GFriction=GFriction+Friction;
  GallFriction=GallFriction+allFriction;
  GVelocity=GVelocity+Velocity;
  GallVelocity=GallVelocity+allVelocity;
end
GDensity=GDensity/length(files);
GallLength=GallLength/length(files);
GCurv=GCurv/length(files);
GFriction=GFriction/length(files);
GallFriction=GallFriction/length(files);
GVelocity=GVelocity/length(files);
GallVelocity=GallVelocity/length(files);
if strcmp(p.Results.File,'off')
    disp('not saving to file holding in memory')
else
    fout=[p.Results.File,'.mat'];
    save(fout,'GDensity','GCurv','GFriction','GallFriction','GVelocity','GallVelocity','GallLength')
end
    
