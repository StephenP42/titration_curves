from math import log10
from numpy import arange
from pylab import *
import pandas as pd

#Исходные данные
print('Al^3+ - титруемое', '\n'
      'Y^4- - титрант')
lgb = 16.5 #Логарифм константы устойчивости
b = float(10**lgb) #Константа устойчивости
ct = 0.1 #Концентрация титранта
co = 0.1 #Концентрация титруемого
ph = 12.0 #Без комментариев
h = float(10**-ph)

listpka = [2, 2.67, 6.16, 10.26]
listka = []
for i in listpka:
    kael = 10**-i
    listka.append(kael)

ka = 1
for j in listka:
    ka = ka*j

moldolTRIL = ka/(h**4+listka[0]*h**3+listka[0]*listka[1]*h**2+listka[0]*listka[1]*listka[2]*h**1+listka[0]*listka[1]*listka[2]*listka[3])
moldolME = 0.9

#вычислим шаг и pme интервал
pMeINThigh = lgb+log10(moldolTRIL)+log10(1.7-1)
pMeINTlow = -log10(co)
step = (-pMeINTlow + pMeINThigh)/1401

#расчёт
pMe = arange(pMeINThigh, pMeINTlow, -step)
me = 10**-pMe
tchislitel = (ct*co/me)-ct-b*moldolTRIL*me*co+b*moldolTRIL*moldolME*ct*co
znamenatel = b*moldolTRIL*moldolME*ct*co+b*moldolME*me*co+co
ft = tchislitel/znamenatel

#таблица
columns = {'pMe': pMe, 'Me': me, 'Titration Factor': ft}
df = pd.DataFrame(columns)
print(df)
#df.to_csv('KO.csv')

figure()
plot(ft, pMe)
xlim(0, 2)
ylim(0, 18)
xlabel('Фактор оттитрованности')
ylabel('pMe')
title('Кривая комплексонометрического титрования')
show()