def sia_approx(timestep):
    return 104.3 - 2.394*timestep + 28*timestep**0.5732

def sia_approx_low(timestep):
    return 44.2 - 3.14*timestep + 22.1*timestep**0.6584

L = 3.554*20

for i in [t+0.5 for t in range(101)]:
    print(sia_approx_low(i)/L**3, end = ', ')