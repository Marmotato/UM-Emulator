 CIR_total1=Final_responsel1+Final_responsen1+Final_responsen11+Final_responses1;
 CIR_total2=Final_responsel2+Final_responsen2+Final_responsen12+Final_responses2;
 CIR_total3=Final_responsel3+Final_responsen3+Final_responsen13+Final_responses3;
 CIR_total4=Final_responsel4+Final_responsen4+Final_responsen14+Final_responses4;
 CIR_total5=Final_responsel5+Final_responsen5+Final_responsen15+Final_responses5;
 
 %Figura que compara el CIR total en los 5 escenarios
 figure
 plot3(time1,D1*ones(size(time1)),(CIR_total1(1:length(time1))));
 grid on
 hold on
 xlabel('Time (s)') 
 ylabel('Distance (m)')
 zlabel('CIR')
 plot3(time1,D2*ones(size(time1)),(CIR_total2(1:length(time1))));
 plot3(time1,D3*ones(size(time1)),(CIR_total3(1:length(time1))));
 plot3(time1,D4*ones(size(time1)),(CIR_total4(1:length(time1))));
 plot3(time1,D5*ones(size(time1)),(CIR_total5(1:length(time1))));
 legend([{'Scenario 1'},{'Scenario 2'},{'Scenario 3'},{'Scenario 4'},{'Scenario 5'}],'Location','east')
 
 
 %Figuras que comparan el nLoS (2 paredes) en los 5 escenarios
 CIR_nlos1=Final_responsen1+Final_responsen11;
 CIR_nlos2=Final_responsen2+Final_responsen12;
 CIR_nlos3=Final_responsen3+Final_responsen13;
 CIR_nlos4=Final_responsen4+Final_responsen14;
 CIR_nlos5=Final_responsen5+Final_responsen15;
 figure
 hold on
 grid on
 xlabel('Time (s)') 
 ylabel('CIR')
 plot(time1, CIR_nlos1(1:length(time1)));
 plot(time1, CIR_nlos2(1:length(time1)));
 plot(time1, CIR_nlos3(1:length(time1)));
 plot(time1, CIR_nlos4(1:length(time1)));
 plot(time1, CIR_nlos5(1:length(time1)));
 legend([{'Scenario 1'},{'Scenario 2'},{'Scenario 3'},{'Scenario 4'},{'Scenario 5'}],'Location','east')
   
 
 %Figuras que comparan el scatter N=40 en los 5 escenarios
 figure
 hold on
 grid on
 xlabel('Time (s)') 
 ylabel('CIR')
 plot(time1, Final_responses1(1:length(time1)));
 plot(time1, Final_responses2(1:length(time1)));
 plot(time1, Final_responses3(1:length(time1)));
 plot(time1, Final_responses4(1:length(time1)));
 plot(time1, Final_responses5(1:length(time1)));
 legend([{'Scenario 1'},{'Scenario 2'},{'Scenario 3'},{'Scenario 4'},{'Scenario 5'}],'Location','east')
      
 
 
 %Figura que compara el CIR ideal con el CIR total minero escenario 1
 CIR_totali=Final_responseli+Final_responseni+Final_responsen1i; 
 figure
 hold on
 grid on
 plot(time1, CIR_totali(1:length(time1)));
 plot(time1, CIR_total1(1:length(time1)));
 legend([{'CIR total typical indoor environment'},{'CIR total scenario 1'}],'Location','east')

 
 %Figura que compara el scatter para valores de N en el escenario 1
 figure
 hold on
 grid on
 plot(time1, Final_responses20(1:length(time1)));
 plot(time1, Final_responses1(1:length(time1)));
 plot(time1, Final_responses80(1:length(time1)));
 plot(time1, Final_responses120(1:length(time1)));
 plot(time1, Final_responses160(1:length(time1)));
 plot(time1, Final_responses200(1:length(time1)));
 plot(time1, Final_responses250(1:length(time1)));
 plot(time1, Final_responses300(1:length(time1)));
 legend([{'N=20'},{'N=40'},{'N=80'},{'N=120'},{'N=160'},{'N=200'},{'N=250'},{'N=300'}],'Location','east')

 %Figura que compara nlos+sca del minero con nlos del ideal
 CIR_nloside=Final_responseni+Final_responsen1i;
 CIR_nlossca=Final_responsen1+Final_responsen11+Final_responses1;
 figure
 figure
 hold on
 grid on
 plot(time1, CIR_nloside(1:length(time1)));
 plot(time1, CIR_nlossca(1:length(time1)));
 legend([{'CIR de la componente no-LoS escenario ideal'},{'CIR de la componente no-LoS+scattering escenario 1'}],'Location','east')
  