from numpy import arange
from pylab import *
from matplotlib import mlab
import pandas as pd

#Титруемое
Eo = 0.17
n1 = 2.0
m1 = 4.0
print("H2SO3 + H20 - 2e = SO42- + 4H+") #Полуреакция окисления

#Титрант
Et = 1.45
n2 = 5.0
m2 = 8.0
print("MnO4- + 8H+ + 5e- = Mn2+ + 4H2O") #Полуреакция восстановления

pH = 2.0

#Используемые константы
f = 96485.33
r = 8.31
T = 298.0


tn = 0.001 #фт начало
tk = 1.999 #фт конец


eN = Eo - ((0.0592*m1*pH)/n1) + log10(tn/(1.0-tn))*0.0592/n1 #потенциал раствора в начале
eK = Et - ((0.0592*m2*pH)/n2) + log10((tk-1)/tk)*0.0592/n2 #потенциал раствора в конце
step = (eK - eN)/1400

#расчёты
E = arange(eN, eK+1, step)
tau = (1/(1+10**(-pH*m1)*exp(((Eo-E)*n1*f)/(r*T))))/(1-(1/(1+10**(-pH*n2)*exp(((Et-E)*n2*f)/(r*T)))))

#таблица
columns = {'E': E, 'Titration Factor': tau}
df = pd.DataFrame(columns)
print(df)

#вывод кривой
figure()
plot(tau, E)
xlim(0, 2)
xlabel('Фактор оттитрованности')
ylabel('E')
title('Кривая окислительно-восстановительного титрования')
show()