import numpy as np
from kivy.app import App
from kivy.lang import Builder
from kivy.properties import NumericProperty
from kivy.uix.screenmanager import Screen, ScreenManager
from matplotlib import pyplot, style
from scipy.optimize import fsolve

style.use('ggplot')

np.seterr(divide='ignore', invalid='ignore')

class MainWindow(Screen):
    F = NumericProperty(0)
    zf = NumericProperty(0)
    alpha = NumericProperty(0)
    xd = NumericProperty(0)
    xb = NumericProperty(0)
    R = NumericProperty(0)
    q = NumericProperty(0)

    def mccabe_thiele(self):
        xd = float(self.ids.xd.text)
        xb = float(self.ids.xb.text)
        a = float(self.ids.alpha.text)
        zf = float(self.ids.zf.text)
        F = float(self.ids.F.text)
        q = float(self.ids.q.text)
        R = float(self.ids.R.text)

        def dbf(f):
            return [xd * f[0] + xb * f[1] - zf * F, f[0] + f[1] - F]

        [D, B] = fsolve(dbf, [30, 20])

        Lr = R * D
        Ls = Lr + q * F
        Vr = Lr + D
        Vs = Vr + (q - 1) * F

        def x_eq(x):
            # x on the equlibrium curve
            return x / (a * (1 - x) + x)

        def rec_opline(x):
            # rectifying section operating line: y = (Lr/Vr)*x + (D*xd/Vr)
            return (Lr / Vr) * x + (D * xd / Vr)

        def strip_opline(x):
            # stripping section operating line: y = (Ls/Vs)*x - (B*xb/Vs)
            return (Ls / Vs) * x - (B * xb / Vs)

        # intersection point of rectifying opline and stripping opline
        def inter_pt(p):
            return [(Lr / Vr) * p[0] + (D * xd / Vr) - p[1], (Ls / Vs) * p[0] - (B * xb / Vs) - p[1]]

        [xq, yq] = fsolve(inter_pt, [0.5, 0.5])

        # y-x equilibrium curve
        x = np.linspace(0, 1, 10000)
        y = a * x / (1 + x * (a - 1))

        pyplot.figure(figsize=(5, 5))
        pyplot.suptitle("McCabe-Thiele Plot")
        pyplot.xlabel('x')
        pyplot.ylabel('y')
        pyplot.plot(x, y)
        pyplot.plot(x, x, '--')
        pyplot.xlim(0, 1)
        pyplot.ylim(0, 1)

        # rectifying section operating line
        x, y = [xq, xd], [yq, xd]
        pyplot.plot(x, y, label='rectifying section', color='b')

        # stripping section operating line
        x, y = [xq, xb], [yq, xb]
        pyplot.plot(x, y, label='stripping section', color='g')
        pyplot.legend(loc='best')

        x0, y0 = xd, xd
        for i in range(1, 100):
            x1, y1 = x_eq(y0), y0
            pyplot.plot([x0, x1], [y0, y1], color='r')
            if x1 > xq:
                x2, y2 = x1, rec_opline(x1)
            if x1 < xq:
                x2, y2 = x1, strip_opline(x1)
            if (x2, y2) < (xb, xb):
                pyplot.plot([x1, x2], [y1, x2], color='r')
            else:
                pyplot.plot([x1, x2], [y1, y2], color='r')
            x0, y0 = x2, y2
            if y2 < x2:
                break
        pyplot.title("Number of stages = infinity", size=10) if i == 99 \
            else pyplot.title("Number of stages = %d" % i, size=10)
        pyplot.show()


class SecondWindow(Screen):
    alpha = NumericProperty(0)
    xd = NumericProperty(0)
    xb = NumericProperty(0)

    def total_reflux(self):
        xd = float(self.ids.xd.text)
        xb = float(self.ids.xb.text)
        a = float(self.ids.alpha.text)

        def x_eq(x):
            # x on the equlibrium curve
            return x / (a * (1 - x) + x)

        # y-x equilibrium curve
        x = np.linspace(0, 1, 10000)
        y = a * x / (1 + x * (a - 1))

        pyplot.figure(figsize=(5, 5))
        pyplot.suptitle("McCabe-Thiele Plot - Total Reflux")
        pyplot.xlabel('x')
        pyplot.ylabel('y')
        pyplot.plot(x, y)
        pyplot.plot(x, x)
        pyplot.xlim(0, 1)
        pyplot.ylim(0, 1)

        x0, y0, xb = xd, xd, xb
        for i in range(1, 100):
            x1, y1 = x_eq(y0), x0
            pyplot.plot([x0, x1], [y0, y1], color='r')
            x2, y2 = x1, x1
            pyplot.plot([x1, x2], [y1, x2], color='r')
            x0, y0 = x2, y2
            if x2 < xb:
                break

        pyplot.title("Number of stages = infinity", size=10) if i == 99 \
            else pyplot.title("Number of stages = %d" % i, size=10)
        pyplot.show()


class WindowManager(ScreenManager):
    pass


kv = Builder.load_file("mccabethiele.kv")


class mccabethiele(App):
    def build(self):
        return kv


mccabethiele().run()
