import numpy as np
import matplotlib.pyplot as plt

x = dict()
x[0] = [9]*2+[4]*41+[3]*111+[2]*296+[1]*1157+[6]*9+[5]*26+[7]*9+[11]*1+[8]*1+[10]*1
x[9] = [24]*1+[10]*10+[9]*19+[8]*24+[7]*30+[6]*50+[5]*39+[4]*52+[3]*75+[2]*104+[1]*257+[14]*1+[11]*9+[13]*8+[12]*8+[19]*1+[16]*3+[15]*3+[25]*1
x[19] = [24]*1+[16]*9+[13]*8+[11]*9+[10]*11+[9]*20+[8]*15+[6]*38+[5]*24+[3]*23+[2]*39+[1]*99+[22]*3+[21]*2+[12]*14+[4]*26+[28]*1+[14]*5+[7]*24+[30]*1+[27]*1+[23]*1+[15]*5+[18]*5+[20]*4+[17]*3+[25]*1+[19]*3+[26]*1+[35]*1
x[29] = [29]*3+[24]*2+[15]*7+[12]*14+[11]*6+[9]*13+[8]*11+[6]*26+[5]*14+[3]*12+[1]*49+[35]*1+[27]*3+[22]*2+[13]*9+[4]*15+[28]*2+[18]*6+[14]*5+[7]*20+[31]*1+[23]*2+[2]*13+[17]*7+[21]*1+[19]*4+[16]*6+[38]*1+[20]*2+[10]*5+[25]*2+[34]*1+[26]*3+[51]*1+[42]*1+[41]*1
x[25] = [24]*1+[16]*8+[15]*8+[13]*9+[12]*14+[10]*6+[9]*15+[8]*9+[6]*31+[5]*18+[3]*16+[1]*65+[35]*1+[27]*3+[22]*2+[4]*14+[28]*2+[14]*4+[11]*9+[7]*25+[31]*1+[23]*2+[2]*21+[17]*6+[19]*6+[32]*1+[18]*5+[29]*1+[20]*1+[25]*2+[34]*1+[26]*2+[42]*1+[41]*1
x[49] = [44]*2+[41]*1+[24]*2+[14]*4+[3]*5+[1]*19+[36]*3+[22]*4+[82]*1+[9]*6+[8]*9+[5]*5+[31]*2+[30]*1+[28]*3+[23]*6+[7]*10+[19]*2+[17]*3+[13]*3+[35]*2+[16]*2+[6]*14+[4]*4+[2]*3+[68]*1+[18]*4+[32]*2+[26]*2+[15]*4+[54]*1+[46]*1+[29]*2+[20]*3+[10]*3+[47]*1+[27]*2+[37]*1+[25]*1+[40]*1+[12]*5+[11]*1+[73]*1+[21]*1+[76]*1+[50]*1+[66]*1
x[75] = [43]*1+[41]*2+[27]*4+[15]*3+[1]*9+[44]*1+[39]*2+[22]*3+[104]*1+[14]*3+[9]*5+[55]*1+[30]*2+[28]*5+[65]*1+[20]*3+[13]*2+[70]*1+[8]*6+[7]*5+[6]*9+[97]*1+[18]*2+[40]*1+[33]*1+[23]*2+[54]*1+[46]*1+[21]*1+[34]*1+[29]*3+[47]*1+[38]*1+[37]*1+[17]*1+[64]*1+[11]*1+[112]*1+[12]*1+[53]*1+[51]*1+[79]*1+[10]*1+[3]*1+[77]*1+[31]*1+[4]*1+[89]*1+[32]*1
x[99] = [44]*1+[41]*3+[27]*2+[15]*3+[65]*3+[39]*1+[22]*1+[1]*9+[104]*1+[14]*2+[9]*4+[69]*1+[57]*2+[28]*1+[20]*3+[13]*1+[106]*1+[6]*7+[97]*1+[18]*1+[40]*1+[23]*2+[8]*3+[75]*3+[46]*2+[37]*3+[34]*2+[29]*2+[7]*5+[47]*1+[113]*1+[42]*1+[53]*1+[51]*1+[119]*1+[31]*1+[19]*1+[88]*1+[32]*1

xj = dict()
hist_data = open("hist_kinetics_data.txt", 'r')
exec(hist_data.read())
hist_data.close()

xp = dict()
hist_data = open("output_100_1.txt", 'r')
exec(hist_data.read())
hist_data.close()

def set_xs(t):
    return x[t], xj[t], 0#xp[t]

def show_hist(x, xj, xp, show, show_bars=0):
    show_x, show_xj, show_xp = show
    # the histogram of the data
    if show_xp:
        n, bins, patches = plt.hist(xp, int(max(xp)), range = (0.5, max(xp)+0.5), facecolor='r', alpha = 0.75)
    if show_x:
        n, bins, patches = plt.hist(x, int(max(x)//2), range = (0.5, max(x)+0.5), facecolor='g', alpha = 0.75)
        print(' '.join(map(str, n)))
    if show_xj:
        plt.plot(range(1, 121), xj, color='b')

    if show_bars:
        ##### from StackOverflow
        bincenters = 0.5*(bins[1:]+bins[:-1])
        menStd     = 2 * n ** 0.5
        width      = 0.05
        plt.bar(bincenters, n, width=width, color='g', alpha = 0, yerr=menStd)

    plt.xlabel('???????????? ????????????, (??????-???? ??????????????)')
    plt.ylabel('??????-???? ?????????????? ???????????? ??????????????')
    # plt.xlim(0, 100)
    # plt.ylim(0, 0.8)
    # plt.grid(True)
    # plt.show()

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

X, Xj, Xp = set_xs(49)
show_hist(X, Xj, Xp, (1, 1, 0))

# print_size_growth(x, xj, xp, (0, 1, 0))