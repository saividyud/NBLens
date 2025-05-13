import matplotlib.pyplot as plt

class LinePlot:
    def __init__(self):
        # Each instance gets its own figure and axes
        print('Hello')
    
    def plot_line(self, x, y, label=None):
        self.fig, self.ax = plt.subplots()
        self.ax.plot(x, y, label=label)
        if label:
            self.ax.legend()

plot1 = LinePlot()
plot2 = LinePlot()

plot1.plot_line([1, 2, 3], [1, 4, 9], label='squares')
plot2.plot_line([1, 2, 3], [1, 2, 3], label='linear')

plot1.ax.scatter([1, 2], [1, 4], color='red', label='scatter')

plt.show()