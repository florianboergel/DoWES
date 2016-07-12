import matplotlib.pyplot as plt
import numpy as np

x = np.arange(1,8760)

def weib(x,n,a):
    return (a / n) * (x / n)**(a - 1) * np.exp(-(x / n)**a)

count, bins, ignored = plt.hist(np.random.weibull(5.,1000))
plt.plot(x, weib(x, 2., 9.))
plt.show()
