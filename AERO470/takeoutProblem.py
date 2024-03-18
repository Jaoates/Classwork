import numpy as np
import matplotlib.pyplot as plt

class Simulation:
    def __init__(self, T):
        # System States
        self.N = 0

        # Simulation Variables
        self.clock = 0

        # Event List
        self.t_arrival = self.generate_arrival()
        self.t_depart = T

        # Statistical Counters
        self.N_arrivals = 0
        self.N_departs = 0
        self.total_wait = 0.0
        self.total_wait_vec = []

    def advance_time(self):
        t_event = min(self.t_arrival, self.t_depart)
        self.total_wait += self.N*(t_event - self.clock)
        self.clock = t_event
        if t_event == self.t_arrival:
            self.handle_arrival()
        else:
            self.handle_departure()

    def handle_arrival(self):
        self.N_arrivals += 1
        self.N += 1
        if self.N <= 1:
            self.t_depart = self.clock + self.generate_service()

        self.t_arrival = self.clock + self.generate_arrival()

    def handle_departure(self):
        self.N_departs += 1
        self.N -= 1
        if self.N <= 0:
            self.t_depart = float("inf")

    def generate_arrival(self):
        return np.random.exponential(1./3)

    def generate_service(self):
        return np.random.exponential(1./4)

np.random.seed(0)
sim = Simulation(float("inf"))

tarr = []
tdep = []
Ncust = []
t = []
for i in range(100):
    sim.advance_time()
    tarr.append(sim.t_arrival)
    tdep.append(sim.t_depart)
    t.append(sim.clock)
    Ncust.append(sim.N)
    sim.total_wait_vec.append(sim.total_wait)

print(f"{sim.N_departs} customers were served")
print(f"The total wait time is {sim.total_wait} minutes")

fig, ax = plt.subplots()
ax.plot(t, sim.total_wait_vec)
plt.show()

fig, ax = plt.subplots()
ax.plot(t, Ncust)
plt.show()
