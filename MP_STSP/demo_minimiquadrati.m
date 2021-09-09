first=1
last=n_cols

a=val_grid_utm_norm(first,2)
b=val_grid_utm_norm(last,2)
x=val_grid_utm_norm(1:last,2)
y=val_grid_utm_norm(1:last,3)

n=1
%valutazionecoeff.dell'approssimanteaiminimiquadrati
%”pn”di”f”digrado”n”
coeff=polyfit(x,y,n);
%valore”pn”nelleascisse”x”
z=polyval(coeff,x);
%errore||f-pn||2
err2=norm(z-y,2);
fprintf('\n\tn:%2.0fnorma2:%1.2e',n,err2);
%graficodelpolinomioaiminimiquadratidigrado”n”
ht=1/10000;u=a:ht:b;
v=polyval(coeff,u);
clf;
plot(x,y,'r.');
hold on;
plot(u,v,'k-','LineWidth',2);
titlestr=strcat('Minimiquadratidigrado:',num2str(n));
title(titlestr);
legend('Dati','ApprossimanteMinimiquadrati');
hold off;

%theta_grid=atan2(max(val_grid_utm_norm(1:last,3)),max(val_grid_utm_norm(1:last,2)))
theta_grid=atan2(val_grid_utm_norm(last,3),val_grid_utm_norm(last,2))