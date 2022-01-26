import numpy as np
import matplotlib.pyplot as plt

x = dict()
x[1] = [9]*2+[4]*41+[3]*111+[2]*296+[1]*1157+[6]*9+[5]*26+[7]*9+[11]*1+[8]*1+[10]*1
x[25] = [24]*1+[16]*8+[15]*8+[13]*9+[12]*14+[10]*6+[9]*15+[8]*9+[6]*31+[5]*18+[3]*16+[1]*65+[35]*1+[27]*3+[22]*2+[4]*14+[28]*2+[14]*4+[11]*9+[7]*25+[31]*1+[23]*2+[2]*21+[17]*6+[19]*6+[32]*1+[18]*5+[29]*1+[20]*1+[25]*2+[34]*1+[26]*2+[42]*1+[41]*1
x[50] = [44]*2+[41]*1+[24]*2+[14]*4+[3]*5+[1]*17+[35]*3+[22]*4+[81]*1+[9]*6+[8]*9+[5]*5+[31]*2+[30]*1+[28]*3+[23]*6+[7]*10+[19]*1+[17]*3+[13]*3+[36]*2+[16]*2+[6]*14+[4]*4+[2]*3+[69]*1+[18]*4+[33]*1+[26]*2+[15]*4+[54]*1+[46]*1+[29]*2+[20]*3+[10]*3+[47]*1+[27]*2+[37]*1+[25]*1+[40]*1+[12]*5+[11]*1+[73]*1+[21]*2+[76]*1+[50]*1+[66]*1+[32]*1
x[75] = [43]*1+[41]*2+[27]*4+[15]*3+[1]*9+[44]*1+[39]*2+[22]*3+[104]*1+[14]*3+[9]*5+[55]*1+[30]*2+[28]*5+[65]*1+[20]*3+[13]*2+[70]*1+[8]*6+[7]*5+[6]*9+[97]*1+[18]*2+[40]*1+[33]*1+[23]*2+[54]*1+[46]*1+[21]*1+[34]*1+[29]*3+[47]*1+[38]*1+[37]*1+[17]*1+[64]*1+[11]*1+[112]*1+[12]*1+[53]*1+[51]*1+[79]*1+[10]*1+[3]*1+[77]*1+[31]*1+[4]*1+[89]*1+[32]*1
x[100] = [44]*1+[41]*3+[27]*2+[15]*3+[105]*1+[22]*1+[104]*1+[14]*2+[9]*4+[69]*1+[58]*1+[65]*2+[28]*1+[20]*3+[13]*1+[1]*4+[106]*1+[6]*7+[97]*1+[18]*1+[40]*1+[23]*2+[8]*3+[75]*2+[46]*2+[37]*3+[34]*2+[29]*2+[7]*5+[47]*1+[113]*1+[42]*1+[53]*1+[51]*1+[57]*1+[119]*1+[77]*1+[31]*1+[19]*1+[89]*1+[32]*1

xj = dict()
hist_data = open("hist_kinetics_data.txt", 'r')
exec(hist_data.read())
hist_data.close()

xp = dict()
hist_data = open("output_100_1.txt", 'r')
exec(hist_data.read())
hist_data.close()

def set_xs(t):
    return x[t], xj[t], xp[t]

def show_hist(x, xj, xp, show, show_bars=0):
    show_x, show_xj, show_xp = show
    # the histogram of the data
    if show_xp:
        n, bins, patches = plt.hist(xp, int(max(xp)), range = (0.5, max(xp)+0.5), facecolor='r', alpha = 0.75)
    if show_x:
        n, bins, patches = plt.hist(x, int(max(x)), range = (0.5, max(x)+0.5), facecolor='g', alpha = 0.75)
    if show_xj:
        plt.plot(range(1, 121), xj, color='b')

    if show_bars:
        ##### from StackOverflow
        bincenters = 0.5*(bins[1:]+bins[:-1])
        menStd     = 2 * n ** 0.5
        width      = 0.05
        plt.bar(bincenters, n, width=width, color='g', alpha = 0, yerr=menStd)

    plt.xlabel('Размер пузыря, (кол-во ксенона)')
    plt.ylabel('кол-во пузырей такого размера')
    # plt.xlim(0, 100)
    # plt.ylim(0, 0.8)
    # plt.grid(True)
    plt.show()

def print_size_growth(x, xj, xp, show):
    show_x, show_xj, show_xp = show
    for t in x.keys():
        x_val = np.mean(x[t])
        xp_val = np.mean(xp[t])
        xj_val = np.sum(xj[t]*np.array(range(1, 121))) / np.sum(xj[t])

        to_print = []
        if show_x:
            to_print.append(x_val)
        if show_xj:
            to_print.append(xj_val)
        if show_xp:
            to_print.append(xp_val)
        to_print = list(map(str, to_print))

        print(t, ', '.join(to_print), sep = ', ')

X, Xj, Xp = set_xs(50)
show_hist(X, Xj, Xp, (1, 1, 0))

# print_size_growth(x, xj, xp, (0, 1, 0))