import numpy as np
import matplotlib.pyplot as plt
# from scipy import linalg as spl
from scipy import integrate as int
import copy as copy

class Event:
    def __init__(self,name,time,duration,value):
        self.name = name
        self.time = time
        self.duration = duration
        self.value = value

    
    def __repr__(self):
        return f"Event:{self.name}"

class Schedule:
    events = [Event("nullEvent",0,0,0)]
    value = 0
    def __repr__(self):
        return f"Schedule with: value: {self.value} events: {[t.name for t in self.events]}"
    def updateVal(self):
        self.value = sum([t.value for t in self.events])
    def getNextFree(self):
        return self.events[-1].time + self.events[-1].duration
    def append(self,event):
        assert(isinstance(event,Event))
        self.events = self.events + [event]
        self.updateVal()


times = np.array([12,3,7,21,23,24,1,17,13,20])
durations = np.array([4,4,2,3,3,1,3,2,1,4])
values = np.array([5,8,8,9,5,7,5,6,6,5])

n = len(times)
events = []
for i in range(n):
    events.append(Event(i,times[i],durations[i],values[i]))
    # print(events[i])

# events.append(Event("nullEvent",0,0,0))

events = sorted(events,key=lambda x: x.time)
# sched = Schedule()
# sched.events = events
# sched.updateVal()
# print(sched)
schedules = [Schedule()]
for e in events:
    sa =[] # schedules to add
    for s in schedules:
        if s.getNextFree() <= e.time:
            sc = copy.deepcopy(s) # schedule copy, somehow this doesn't deepcopy the actual events list tho :/
            # print()
            # print(s)
            sc.append(e)
            # print(s)
            # print(sc)
            # hmm = id(sc) == id(s)
            sa.append(sc)
    schedules.extend(sa)

schedules = sorted(schedules,key = lambda x: x.value)
for s in schedules: print(s)

print()
print(f"There are {len(schedules)} feasible schedules")
print(f"The best schedule is - {schedules[-1]}")
