import numpy as np

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

    def advance_time(self):
        t_event = min(self.t_arrival, self.t_depart)
        pass

    def handle_arrival(self):
        pass

    def handle_departure(self):
        pass

    def generate_arrival(self):
        return np.random.exponential(1./3)

    def generate_service(self):
        return np.random.exponential(1./4)

np.random.seed(0)
s = Simulation(float('inf'))

for i in range(100):
    s.advance_time()