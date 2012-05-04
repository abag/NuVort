%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_plot.m to create snapshots of filament
function vortex_anim(filenumbers)
eps=0;
figure('visible','off');
for i=filenumbers
  vortex_plot(i)
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
  
