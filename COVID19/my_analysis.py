
import numpy
from srxraylib.plot.gol import set_qt, plot
set_qt()



urlDeaths = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv"
urlCases = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv"



a = numpy.genfromtxt(urlCases,delimiter=' ,',dtype=str)
b = numpy.genfromtxt(urlDeaths,delimiter=' ,',dtype=str)

print(a[0])
print(a[19])
print(a[17])

dates=a[0].split(",")

italy=a[17].split(",")
italyD=b[17].split(",")

spain=a[19].split(",")
spainD=b[19].split(",")

california=a[101].split(",")
californiaD=b[101].split(",")


uptoday=numpy.linspace(-len(spain),0,len(spain)+1)
# print(len(dates), len(spain), len(italy))
for i in range(21):
    print(dates[-i], spain[-i], italy[-i], california[-i])
    # print("s   %s   %s "%(titles[i],spain[i],italy[i]))


# import matplotlib.pylab as plt
# plt.xlabel('Date')
# plt.ylabel('XX')
# plt.bar(uptoday[-29:],numpy.array(italy[-29:],dtype=float),width=0.25)
# plt.bar(uptoday[-29:]-9+0.3,numpy.array(spain[-29:],dtype=float),width=0.25)
# plt.show()

shift=7
t = uptoday[-29:]
s = numpy.array(spain[-29:],dtype=float).copy()
sD = numpy.array(spainD[-29:],dtype=float).copy()

i = numpy.array(italy[-29:],dtype=float).copy()
iD = numpy.array(italyD[-29:],dtype=float).copy()

us = numpy.array(italy[-29:],dtype=float).copy()
usD = numpy.array(italyD[-29:],dtype=float).copy()

plot(t,i,
     t-shift,s,
     t,iD,
     t-shift,sD,
     xtitle="Date up today",legend=["Italy Cases","Spain Cases","Italy Deaths","Spain Deaths"],
     title="Spain = Italy - %d"%shift,ylog=1,marker=['o','o','x','x'],show=0)



t = t
y4 = 2**(t/4)
y4 = y4 / y4[-1]

y2 = 2**(t/3)
y2 = y2 / y2[-1]

plot(
     t,s,
     t,sD,
     t,y4 * s[-1],
     t,y2 * sD[-1],
     xtitle="Date up today",legend=["Spain Cases","Spain Deaths","",""],
     title="Spain",ylog=1)

t = numpy.concatenate((t,numpy.array([1,2])))
y4 = 2**(t/4)
y4 = y4 / y4[-3]

y2 = 2**(t/3)
y2 = y2 / y2[-3]

print("Spain 22/3/2020 %d cases %d deaths"%(s[-1],sD[-1]))
# print("Spain fit %d cases %d deaths"%(s[-1],sD[-1]))
print("Spain prediction %d cases %d deaths"%((y4 * s[-1])[-2],(y2 * sD[-1])[-2]))
print("Spain prediction %d cases %d deaths"%((y4 * s[-1])[-1],(y2 * sD[-1])[-1]))


# Spain 22/3/2020 28768 cases 1772 deaths
# Spain prediction 34211 cases 2232 deaths
# Spain prediction 40684 cases 2812 deaths