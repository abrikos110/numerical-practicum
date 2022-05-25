import numpy
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def runge_kutta_2(x0, y0, f, h=0.01, n=100):
    # y(x0) = y0; y'(x) = f(x, y)
    y = [y0]
    x = [x0]
    for i in range(n):
        k1 = h * f(x[-1], y[-1])
        k2 = h * f(x[-1] + h, y[-1] + k1)
        yi = y[-1] + (k1 + k2) / 2
        xi = x[-1] + h
        y.append(yi)
        x.append(xi)
    return numpy.array(x), numpy.array(y)

def runge_rule(x0, y0, f, h=0.01, n=100):
    x, y1 = runge_kutta_2(x0, y0.copy(), f, h, n)
    x, y2 = runge_kutta_2(x0, y0.copy(), f, h/2, n*2)
    y2 = y2[::2]
    return abs(y1 - y2).mean()


def update(var):
    def f(val):
        globals()[var] = val
        calc()
    return f

def calc():
    global lines, E2, X
    X = int(X)
    #E2 = mu * E1/E3 - 0.01
    y = numpy.array([E1/E3 * mu, E3/mu, -E2 + mu * E1/E3])
    y += eps
    x = 0
    x, y = runge_kutta_2(x, y, f, h=1/int(nh), n=X*int(nh))
    m = [0, 1]
    for i in range(y.shape[1]):
        lines[i].set_xdata(x)
        lines[i].set_ydata(y[:, i])
        m[0] = min(m[0], y[:, i].min())
        m[1] = max(m[1], y[:, i].max())
    d = abs(m[1] - m[0]) * 0.03
    m[0] -= d
    m[1] += d
    lines[0].axes.set_ylim(*m)
    lines[0].axes.set_xlim(x.min(), x.max())
    plt.gcf().canvas.draw_idle()

if __name__ == '__main__':
    mu = 24.7
    E1 = 0.14
    E2 = 0.185
    E3 = 0.156
    eps = 0.264
    y = numpy.array([E1/E3 * mu, E3/mu, -E2 + mu * E1/E3])
    X = 100
    y += eps
    x = 0
    f = lambda x, y: numpy.array([E1 - y[0]*y[1],
                                  y[1] * (-E2 + y[0] - y[2]),
                                  y[2] * (-E3 + mu*y[1])])
    #print(runge_rule(x, y.copy(), f, h=1/100, n=100)
    #      / runge_rule(x, y.copy(), f, h=1/200, n=200))
    X = 100
    nh = 100
    x, y = runge_kutta_2(x, y, f, h=1/int(nh), n=X*int(nh))

    fig, ax = plt.subplots()
    lines = [plt.plot(x, y[:, i], c=['red', 'green', 'blue'][i])[0] for i in range(y.shape[1])]

    vm = {'eps':1, 'mu':100, 'E1':1, 'E2':1, 'E3':1, 'nh':1000, 'X':1000}
    sliders = [Slider(ax=plt.axes([0.1,0.983-i*0.017,0.8,0.02]), label=v, valmin=0, valmax=vm[v], valinit=eval(v))
               for i, v in enumerate(vm)]
    for sl, v in zip(sliders, vm):
        sl.on_changed(update(v))

plt.show()
