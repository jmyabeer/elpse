load TPD_time_xl_intensity_tpd_Ln_Te.dat
tpd=TPD_time_xl_intensity_tpd_Ln_Te;
nm=length(tpd(:,1));
k1=0;
for k=1:nm
    if tpd(k,3)>0.01
        k1=k1+1;
        tpd1(k1,1)=tpd(k,1);
        tpd1(k1,3)=tpd(k,3);
        tpd1(k1,4)=tpd(k,4);
        etaSRS(k1)=tpd(k,3)/1.e14*(tpd(k,5))^(4./3.)/2377./1;%2.5;
    end
end
    
figure(80)
delete(80)
figure(80)
set(gcf,'paperunits','centimeters')
set(gcf,'paperposition',[6 10 15 12])
set(gcf,'units','centimeters')
set(gcf,'position',[6 10 15 12])

subplot(2,2,1)
plot(tpd(1:nm,1),tpd(1:nm,5))
%axis([0 10 0 1000])
set(gca,'FontSize',10,'XTick',([0 1 2 3 4 5 6 7 8 9 10]))
ylabel('Scale length ({\mu}m)','Fontsize',10) 
xlabel('Time (ns)','Fontsize',10) 

subplot(2,2,2)
plot(tpd(1:nm,1),tpd(1:nm,6))
%axis([0 10 0 6])
set(gca,'FontSize',10,'XTick',([0 1 2 3 4 5 6 7 8 9 10]))
ylabel('Electron temperature (keV)','Fontsize',10) 
xlabel('Time (ns)','Fontsize',10) 

subplot(2,2,3)
plot(tpd1(1:k1,1),tpd1(1:k1,3))
%axis([0 10 0 2.e15])
set(gca,'FontSize',10,'XTick',([0 1 2 3 4 5 6 7 8 9 10]))
ylabel('Intensity (W/cm^2)','Fontsize',10) 
xlabel('Time (ns)','Fontsize',10) 

subplot(2,2,4)
plot(tpd1(1:k1,1),tpd1(1:k1,4),'k')
hold on
plot(tpd1(1:k1,1),etaSRS(1:k1))
hold off
%axis([0 10 0 30])
set(gca,'FontSize',10,'XTick',([0 1 2 3 4 5 6 7 8 9 10]))
ylabel('\eta_{SRS} (blue), \eta_{TPD} (black)','Fontsize',10) 
xlabel('Time (ns)','Fontsize',10) 

%print -dtiff -r100 DiagFig.tif
print -dtiff -r100 EP.tif