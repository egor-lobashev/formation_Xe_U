def sia_approx(timestep):
    return 104.3 - 2.394*timestep + 28*timestep**0.5732

L = 3.556*20

for i in [t+0.5 for t in range(101)]:
    print(sia_approx(i)/L**3, end = ', ')