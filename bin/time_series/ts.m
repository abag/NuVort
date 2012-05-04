%read in the ts file and plot various dignostic information
%if given the option print will print to .eps file rather than screen
function ts(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
  disp('will not print to screen but instead to .eps files')
case 'empty'
  otherwise
  disp('incorrect option, aborting script and printing help:')
  help ts
  return
end
A=load('./data/ts.log');
dims=load('./data/dims.log');
t=A(:,2) ; pcount=A(:,3) ; rcount=A(:,4) ; sep=A(:,5) ; l=A(:,6) ; 
maxu=A(:,7) ; maxdu=A(:,8) ; curv=A(:,9) ; rmcount=A(:,10) ;
density=l/(dims(2)^3) ; spacing=1./sqrt(density);

switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'filament information I')      
end
  subplot(2,2,1)
    plot(t,pcount,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('pcount','FontSize',14)
  subplot(2,2,2)
    plot(t,rcount,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('recon count','FontSize',14)
  subplot(2,2,3)
    plot(t,sep,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('avg sep','FontSize',14)
  subplot(2,2,4)
    plot(t,l,'-g','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('Length (\Lambda)','FontSize',14)
if option=='print'
    disp('printing to filament_information.eps')
    print('-depsc','./filament_information.eps')
end
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'velocity information')      
end
  subplot(2,1,1)
    plot(t,maxu,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max(u)','FontSize',14)
  subplot(2,1,2)
    plot(t,maxdu,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max(du)','FontSize',14)
if option=='print'
  disp('printing to velocity_information.eps')
  print('-depsc','./velocity_information.eps')
end
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'filament information II')      
end
  subplot(2,2,1)
    plot(t,curv,'LineWidth',2,'Color',rgb('Chocolate'));
    set(gca,'FontSize',14);
    xlabel('t','FontSize',14);
    ylabel('curv','FontSize',14);
  subplot(2,2,2)
    plot(t,rmcount,'-','LineWidth',2,'Color',rgb('Indigo'));
    set(gca,'FontSize',14);
    xlabel('t','FontSize',14);
    ylabel('rm count','FontSize',14);
  subplot(2,2,3)
    plot(t,density,'-','LineWidth',2,'Color',rgb('HotPink'));
    set(gca,'FontSize',14);
    xlabel('t','FontSize',14);
    ylabel('line density (L)','FontSize',14);
  subplot(2,2,4)
    plot(t,spacing,'-','LineWidth',2,'Color',rgb('DarkSeaGreen'));
    set(gca,'FontSize',14);
    xlabel('t','FontSize',14);
    ylabel('inter vortex spacing','FontSize',12);
if option=='print'
  disp('printing to filament_information2.eps')
  print('-depsc','./filament_information2.eps')
end
 
