function density_spectrum(filenumbers)
global D t 
A=load('./data/ts.log');
counter=0;
for i=filenumbers
  counter=counter+1;
  t(counter)=A(i,2);
  D(counter)=continuum_box(i,'D',1.5);
end
plot(t,D)

