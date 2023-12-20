%   DEFINICIÓN CÓDIGO

L=6000;%m Muestras tomadas del registro
fs=2048;% Frecuencia de muestreo/sample rate
MyCell= readcell('Trial1.csv'); %Lectura del archivo de los datos de EEG
X=cell2mat(MyCell(:,3)); %Selección del canal2 del EEG, columna 3 
t=(1:1/2048:(553-1/2048));%Vector de tiempo
%Creación de la matriz en donde se realizará el análisis de datos
n=round(length(X)/L);% Redondeo a un número entero, división del total de
%datos del EEG en la columna 3 en segmentaciones de 6000
B=zeros(L,n);% Creación de un vector de ceros de n longitud, 188 columnas
for i=1:n % Bucle para crear una matriz con las segmentaciones de los datos
%del EEG en n columnas
B(:,i)=X(L*(i-1)+1:L*i);
end


%   ACTIVIDAD 1

z=input("Escriba la columna ")
for vectorC=X(z*L-L+1:z*L) 
figure(1)
%Espectro Potencia c188
TvectorC=transpose(vectorC);%trasponer vector
t5=transpose(t:L);%Segmentar el vector tiempo
yvectorC=fft(vectorC);%transformada rápida de Fourier
P2 = abs(yvectorC/L);%valor absoluto (módulo) al    ser una variable compleja
P1 = P2(1:L/2+1);%duplicamos el véctor para que tengan la misma longitud
semilogx(0:(fs/L):(fs/2-fs/L),P1(1:L/2))%ploteamos en una escala semilogaritmica;para obtener el espectro de potencia
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de Potencia, ' + "Columna " + z)
grid on
%Comando para guardar automáticamente las gráficas
grafEP = gcf;
exportgraphics(grafEP,'Gráfica Espectro de Potencia.jpg')
figure(2)
%Gráfica del vector
plot(t5,vectorC)
xlabel('Tiempo (s)')
ylabel('Muestras')
title('Segmentación, ' + "Columna " + z)
grafvectorC = gcf;
exportgraphics(grafvectorC,'Gráfica del vector.jpg')
end

%   ACTIVIDAD 2


for vectorC=X(z*L-L+1:z*L)
figure(3)
VectorSinDC=vectorC-mean(vectorC);
%Espectro Potencia
TVectorSinDC=transpose(VectorSinDC);%trasponer vector
t5=transpose(t:L);%Segmentar el vector tiempo
yVectorSinDC=fft(VectorSinDC);%transformada rápida de Fourier
P2 = abs(yVectorSinDC/L);%valor absoluto (módulo) al    ser una variable compleja
P1 = P2(1:L/2+1);%duplicamos el véctor para que tengan la misma longitud
semilogx(0:(fs/L):(fs/2-fs/L),P1(1:L/2))%ploteamos en una escala semilogaritmica;para obtener el espectro de potencia
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de Potencia sin DC, ' + "Columna " + z)
grid on
grafsinDC = gcf;
exportgraphics(grafsinDC,'Gráfica Espectro de Potencia sin DC.jpg')
figure(4)
%Cálculo de la media 
media=mean(vectorC)
plot(t5, VectorSinDC)
title('Segmentación sin DC, ' + "Columna " + z)
xlabel('Tiempo (s)')
ylabel('Vector Sin DC')
grafsinDC = gcf;
exportgraphics(grafsinDC,'Gráfica del vector sin DC.jpg')
%Desviación Estándar sobre la segmentación
desv_Segmentacion=std(vectorC)
%Desviación Estándar eliminando el nivel DC de la segmentación = Perturbación
desv_SinDC=std(VectorSinDC)
end

for vectorC=X(z*L-L+1:z*L)

% Delta (f <4Hz).

for xfiltradaDelta=filter(HdDelta,vectorC)
figure(5)
% Espectro Potencia.
TDelta=transpose(xfiltradaDelta); % Se traspone el vector.
tDelta=transpose(t:L); % Segmentación del vector tiempo.
yDelta=fft(xfiltradaDelta); % Transformada rápida de Fourier.
P2Delta = abs(yDelta/L); % Valor absoluto (módulo) al ser una variable compleja.
P1Delta = P2Delta(1:L/2+1); % Duplicamos el vector para que tengan la misma longitud.
semilogx(0:(fs/L):(fs/2-fs/L),P1Delta(1:L/2)) % Ploteamos en una escala semilogarítmica, para obtener el espectro de potencia.
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de potencia de la Señal Delta, ' + "Columna " + z)
grid on
grafDelta = gcf;
exportgraphics(grafDelta,'Espectro de potencia de la  Señal Delta.jpg')
figure(6)
% Filtro en función del tiempo.
plot(tDelta,xfiltradaDelta)
xlabel('Tiempo de Delta (s)')
ylabel('Filtro Delta')
title('Señal Delta en función del tiempo, ' + "Columna " + z)
grid on
grafDelta1 = gcf;
exportgraphics(grafDelta1,'Gráfica Señal Delta Tiempo.jpg')
end

% Theta (f 4-7Hz).

for xfiltradaTheta=filter(HdTheta,vectorC)
figure(7)
% Espectro Potencia.
TTheta=transpose(xfiltradaTheta); % Se traspone el vector.
tTheta=transpose(t:L); % Segmentación del vector tiempo.
yTheta=fft(xfiltradaTheta); % Transformada rápida de Fourier.
P2Theta = abs(yTheta/L); % Valor absoluto (módulo) al ser una variable compleja.
P1Theta = P2Theta(1:L/2+1); % Duplicamos el vector para que tengan la misma longitud.
semilogx(0:(fs/L):(fs/2-fs/L),P1Theta(1:L/2)) % Ploteamos en una escala semilogarítmica, para obtener el espectro de potencia.
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de potencia de la Señal Theta, ' + "Columna " + z)
grid on
grafTheta = gcf;
exportgraphics(grafTheta,'Gráfica Espectro de potencia de la Señal Theta.jpg')
figure(8)
% Filtro en función del tiempo.
plot(tTheta,xfiltradaTheta)
xlabel('Tiempo de Theta (s)')
ylabel('Filtro Theta')
title('Señal Theta en función del tiempo, ' + "Columna " + z)
grid on
grafTheta1 = gcf;
exportgraphics(grafTheta1,'Gráfica Señal Theta Tiempo.jpg')
end

% Alfa (f 8-13Hz).

for xfiltradaAlfa=filter(HdAlfa,vectorC)
figure(9)
% Espectro Potencia.
TAlfa=transpose(xfiltradaAlfa); % Se traspone el vector.
tAlfa=transpose(t:L); % Segmentación del vector tiempo.
yAlfa=fft(xfiltradaAlfa); % Transformada rápida de Fourier.
P2Alfa = abs(yAlfa/L); % Valor absoluto (módulo) al ser una variable compleja.
P1Alfa = P2Alfa(1:L/2+1); % Duplicamos el vector para que tengan la misma longitud.
semilogx(0:(fs/L):(fs/2-fs/L),P1Alfa(1:L/2)) % Ploteamos en una escala semilogarítmica, para obtener el espectro de potencia.
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de potencia de la Señal Alfa, ' + "Columna " + z)
grid on
grafAlfa = gcf;
exportgraphics(grafAlfa,'Gráfica Espectro de potencia de la Señal Alfa.jpg')
figure(10)
%Filtro en función del tiempo.
plot(tAlfa,xfiltradaAlfa)
xlabel('Tiempo de Alfa (s)')
ylabel('Filtro Alfa')
title('Señal Alfa en función del tiempo, ' + "Columna " + z)
grid on
grafAlfa1 = gcf;
exportgraphics(grafAlfa1,'Gráfica Señal Alfa Tiempo.jpg')
end

% Beta (f 14-30Hz).

for xfiltradaBeta=filter(HdBeta,vectorC)
figure(11)
% Espectro Potencia.
TBeta=transpose(xfiltradaBeta); % Se traspone el vector.
tBeta=transpose(t:L); % Segmentación del vector tiempo.
yBeta=fft(xfiltradaBeta); % Transformada rápida de Fourier.
P2Beta = abs(yBeta/L); % Valor absoluto (módulo) al ser una variable compleja.
P1Beta = P2Beta(1:L/2+1); % Duplicamos el vector para que tengan la misma longitud.
semilogx(0:(fs/L):(fs/2-fs/L),P1Beta(1:L/2)) % Ploteamos en una escala semilogarítmica, para obtener el espectro de potencia.
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de potencia de la Señal Beta, ' + "Columna " + z)
grid on
grafBeta = gcf;
exportgraphics(grafBeta,'Gráfica Espectro de potencia de la Señal Beta.jpg')
figure(12)
% Filtro en función del tiempo.
plot(tBeta,xfiltradaBeta)
xlabel('Tiempo de Beta (s)')
ylabel('Filtro Beta')
title('Señal Beta en función del tiempo, ' + "Columna " + z)
grid on
grafBeta1 = gcf;
exportgraphics(grafBeta1,'Gráfica Señal Beta Tiempo.jpg')
end

% Gamma (f >30Hz).

for xfiltradaGamma=filter(HdBeta,vectorC)
figure(13)
% Espectro Potencia.
TGamma=transpose(xfiltradaGamma); % Se traspone el vector.
tGamma=transpose(t:L); % Segmentación del vector tiempo.
yGamma=fft(xfiltradaGamma); % Transformada rápida de Fourier.
P2Gamma = abs(yGamma/L); % Valor absoluto (módulo) al ser una variable compleja.
P1Gamma = P2Gamma(1:L/2+1); % Duplicamos el vector para que tengan la misma longitud.
semilogx(0:(fs/L):(fs/2-fs/L),P1Gamma(1:L/2)) % Ploteamos en una escala semilogarítmica, para obtener el espectro de potencia.
xlabel('Frecuencia (Hz)')
ylabel('Magnitud |P1(f)|')
title('Espectro de potencia de la Señal Gamma, ' + "Columna " + z)
grid on
grafGamma = gcf;
exportgraphics(grafGamma,'Gráfica Espectro de potencia de la Señal Gamma.jpg')
figure(14)
% Filtro en función del tiempo.
plot(tGamma,xfiltradaGamma)
xlabel('Tiempo de Gamma (s)')
ylabel('Filtro Gamma')
title('Señal Gamma en función del tiempo, ' + "Columna " + z)
grid on
grafGamma1 = gcf;
exportgraphics(grafGamma1,'Gráfica Señal Gamma Tiempo.jpg')
end

end
