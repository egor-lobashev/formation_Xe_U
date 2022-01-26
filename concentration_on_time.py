f = open('C.txt', 'r')
C = eval(f.readline())
f.close()

def print_concentration_on_time(sizes_to_print):
    for i in range(100):
        try:
            list_to_print = [str(C[(n, 1)]) for n in sizes_to_print]
            print(i, ', '.join(list_to_print), sep = ', ')
        except:
            pass

def C_start():
    C_0 = []
    for n in range(120):
        if (0, n) in C.keys():
            C_0.append(C[(0, n)])
        else:
            C_0.append(0)
    return C_0

# print_concentration_on_time((1))
# print_concentration_on_time((2, 3, 4))

C_0 = C_start()

for n in range(1, len(C_0)):
    print(n, C_0[n], sep = ', ')