from math import log10
from numpy import *
from pylab import *
import pandas as pd

#int
print('2NaBr + Hg(NO3) = Hg2Br2 + NaNO3')
pksp = 22.29
ks = 10**-pksp
co = 0.1
ct = 0.1
x = 2
pHg1 = -log10(ct) # При ФТ больше 1
pHg2 = pksp-1 # При ФТ меньше 1
hg1 = 10**-pHg1
hg2 = 10**-pHg2

#расчёты
step = (pHg2 - pHg1)/1401
phg = arange(pHg1, pHg2, step)
hg = 10**-phg
ft = x*(-x*hg*ct-co*ct+((ks/hg)**(1/x))*ct)/(-co*ct*x-((ks/hg)**(1/x))*co+x*hg*co)

#таблицы
columns = {'Hg+': hg, 'pHg': phg, 'Titration Factor': ft}
df = pd.DataFrame(columns)
print(df)
#df.to_csv('KO.csv')

figure()
plot(ft, phg)
xlabel('Фактор оттитрованности')
xlim(0,2)
ylim(0,26)
ylabel('pHg')
show()