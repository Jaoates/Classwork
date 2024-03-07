import matplotlib.pyplot as plt
import numpy as np

data = {'a': np.arange(50),
        'c': np.random.randint(0, 50, 50),
        'd': np.random.randn(50)}
data['b'] = data['a'] + 10 * np.random.randn(50)
data['d'] = np.abs(data['d']) * 100
data['c'] = data['a']+np.random.randint(0,10)

plt.scatter('a', 'b', c='c', s='d', data=data)
plt.xlabel('entry a')
plt.ylabel('entry b')
plt.show()

theta = np.arange(0,9*np.pi,.1)
r = theta/np.pi
x = np.cos(theta)*r
y = np.sin(theta)*r

plt.scatter(x,y,c=theta,s=(r[-1]-r)*100)
plt.show()

