from numpy import arange
from pylab import *
from matplotlib import mlab
import pandas as pd

#int
c0 = float(0.1)
ct = float(0.1)
pKa = 10
ka = float(10**(-pKa))
x = 1 #при 1 - титруем основание, при 0 - титруем кислоту
pKw = 14.0
kw = float(10**(-pKw))

#расчёты
pHintBeg = 14*x-1**(1-x)*log10(c0)
pHintEnd = (1-x)*14-1**(x)*log10(ct)
step = ( (pHintEnd-pHintBeg)/1401 )
pH = arange(pHintBeg, pHintEnd, step)
h = 10.0**(-pH)
chislitel = ((c0*ct*ka)/(h+ka))+(ct*kw/h)-(c0*ct*x)-(h*ct)
znamenatel = (h*c0)+(c0*ct*(float(1)-float(2)*x))-(c0*kw/h)
ft = chislitel/znamenatel

#таблица
columns = {'pH': pH, 'H+': h, 'Titration Factor': ft}
df = pd.DataFrame(columns)
print(df)
#df.to_csv('KO.csv')

figure()
plot(ft, pH)
xlim(0, 2)
ylim(0, 14)
xlabel('Фактор оттитрованности')
ylabel('pH')
title('Кривая кислотно-основного титрования')
show()



