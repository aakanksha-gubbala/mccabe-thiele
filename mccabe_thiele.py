from tkinter import *

import numpy as np
from matplotlib import pyplot, style
from scipy.optimize import fsolve

style.use('seaborn-white')

root = Tk()
root.title('McCabe Thiele Graphs')
root.geometry('650x300')
root.maxsize(width=650, height=300)


def open():
    top = Toplevel()
    top.title('Total Reflux')
    top.geometry('380x250')

    alpha_label = Label(top, text='alpha', font=('bold, 14'), pady=20)
    alpha_label.grid(row=0, column=0, sticky=W)
    alpha_entry = Entry(top)
    alpha_entry.grid(row=0, column=1)

    xd_label = Label(top, text='Distillate concentration',
                     font=('bold, 14'), pady=20)
    xd_label.grid(row=1, column=0, sticky=W)
    xd_entry = Entry(top)
    xd_entry.grid(row=1, column=1)

    xb_label = Label(top, text='Bottoms concentration',
                     font=('bold, 14'), pady=20)
    xb_label.grid(row=2, column=0, sticky=W)
    xb_entry = Entry(top)
    xb_entry.grid(row=2, column=1)

    def total_reflux():
        xd = float(xd_entry.get())
        xb = float(xb_entry.get())
        a = float(alpha_entry.get())

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
            x1, y1 = x_eq(x0), x0
            pyplot.plot([x0, x1], [y0, y1], color='r')
            x2, y2 = x1, x1
            pyplot.plot([x1, x2], [y1, x2], color='r')
            x0, y0 = x2, y2
            if x2 < xb:
                break

        pyplot.title("Number of stages = infinity", size=10) if i == 99 \
            else pyplot.title("Number of stages = %d" % i, size=10)
        pyplot.show()

    button_run = Button(top, text='Run', font=('bold, 14'),
                        width=10, command=total_reflux)
    button_run.grid(row=4, column=0)

    button_quit = Button(top, text='Close', font=(
        'bold, 14'), width=10, command=top.destroy)
    button_quit.grid(row=4, column=1)


feed_label = Label(root, text='Feed flow rate', font=('bold, 14'), pady=20)
feed_label.grid(row=0, column=0, sticky=W)
feed_entry = Entry(root)
feed_entry.grid(row=0, column=1)

zf_label = Label(root, text='Feed concentration', font=('bold, 14'), pady=20)
zf_label.grid(row=0, column=2, sticky=W)
zf_entry = Entry(root)
zf_entry.grid(row=0, column=3)

alpha_label = Label(root, text='alpha', font=('bold, 14'), pady=20)
alpha_label.grid(row=1, column=0, sticky=W)
alpha_entry = Entry(root)
alpha_entry.grid(row=1, column=1)

xd_label = Label(root, text='Distillate concentration',
                 font=('bold, 14'), pady=20)
xd_label.grid(row=1, column=2, sticky=W)
xd_entry = Entry(root)
xd_entry.grid(row=1, column=3)

r_label = Label(root, text='Reflux Ratio', font=('bold, 14'), pady=20)
r_label.grid(row=2, column=0, sticky=W)
r_entry = Entry(root)
r_entry.grid(row=2, column=1)

xb_label = Label(root, text='Bottoms concentration',
                 font=('bold, 14'), pady=20)
xb_label.grid(row=2, column=2, sticky=W)
xb_entry = Entry(root)
xb_entry.grid(row=2, column=3)

q_label = Label(root, text='Thermal ratio', font=('bold, 14'), pady=20)
q_label.grid(row=3, column=0, sticky=W)
q_entry = Entry(root)
q_entry.grid(row=3, column=1)


def mccabe_thiele():
    xd = float(xd_entry.get())
    xb = float(xb_entry.get())
    zf = float(zf_entry.get())
    F = float(feed_entry.get())
    a = float(alpha_entry.get())
    q = float(q_entry.get())
    R = float(r_entry.get())

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


button_run = Button(root, text='Run', font=('bold, 14'),
                    width=10, command=mccabe_thiele)
button_run.grid(row=4, column=0)

button_tr = Button(root, text='Total Reflux', width=10, command=open)
button_tr.grid(row=4, column=1)

button_quit = Button(root, text='Quit', font=(
    'bold, 14'), width=10, command=root.quit)
button_quit.grid(row=4, column=2)

root.mainloop()
