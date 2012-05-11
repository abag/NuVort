%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_plot.m to create snapshots of filament
function vortex_anim(filenumbers)
eps=0;
figure('visible','off');
for i=filenumbers
  subplot(2,1,1)
    vortex_plot(i,'OverHead','on','LineWidth',2)
  subplot(2,1,2)
    vortex_plot(i,'LineWidth',2)
  if eps==1 
    fOUT=sprintf('data/var%04d.eps',i)
    print('-depsc',fOUT)
  else
    fOUT=sprintf('data/var%04d.png',i)
    print('-dpng',fOUT)
  end 
  %close all open files
  fclose('all');
end
figure('visible','on');
close all
  
