import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from srxraylib.plot.gol import set_qt
set_qt()


def xsin(x,xc):
    g = np.exp(- x ** 2 / 2 / 40 ** 2)
    a = g * np.sin((x - xc) / 8)
    return 0.3 * g + a**2

def animate(i):
    line.set_ydata(F[i, :])


#
x = np.linspace(-100,100,500)
#
#
XC = np.array([-100,-60,-20,20,60,100])


fig, ax = plt.subplots(figsize=(5, 3))

ax.set(xlim=(-100, 100), ylim=(0, 1.5))

t = XC # np.linspace(1, 25, 30)
X2, T2 = np.meshgrid(x, t)

sinT2 = np.sin(2 * np.pi * T2 / T2.max())
sinT2 = xsin(X2,T2)
F = sinT2 # 0.9 * sinT2 * np.sinc(X2 * (1 + sinT2))

line = ax.plot(x, F[0, :], color='k', lw=2)[0]


anim = FuncAnimation(
    fig, animate, interval=100, frames=len(t) - 1)

plt.draw()
plt.show()


anim.save('anholonomy.gif', writer='imagemagick')
# anim.save('anholonomy.mp4')