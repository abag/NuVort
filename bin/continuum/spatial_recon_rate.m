function spatial_recon_rate(dir)
N=32;
load ./data/dims.log
A=load('./data/recon_location.log');
spat_recon_rate(1:N)=0;
if strcmp(dir,'y')     
    disp('recon rate in y direction')
    for i=1:N+1
        xmesh(i)=(i-1)/N*dims(3)-dims(3)/2;
    end
    for j=1000:length(A)
       for i=1:N
           if A(j,3)>xmesh(i) && A(j,3)<xmesh(i+1)
           spat_recon_rate(i)=spat_recon_rate(i)+1;
           end
       end
    end
elseif strcmp(dir,'z')
    disp('recon rate in z direction')
    for i=1:N+1
        xmesh(i)=(i-1)/N*dims(4)-dims(4)/2
    end
    for j=1000:length(A)
       for i=1:N
           if A(j,4)>xmesh(i) && A(j,4)<xmesh(i+1)
           spat_recon_rate(i)=spat_recon_rate(i)+1;
           end
       end
    end
else
    disp('assuming direction is x')
    for i=1:N
        xmesh(i)=(i-1)/N*dims(2)-dims(2)/2
    end
end
xmesh(N+1)=[];
plot(xmesh,spat_recon_rate)