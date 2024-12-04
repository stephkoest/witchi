import json
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import gnuplotlib as gp

class ChiSquarePlotter:
    def __init__(self, json_file, png_file):
        self.json_file = json_file
        self.data = self.load_data()
        self.png_file = png_file

    def load_data(self):
        with open(self.json_file, 'r') as file:
            return json.load(file)

    def plot_distributions_ascii(self):
        labels = []
        distributions = []
        for label, values in self.data.items():
            labels.append(label)
            distributions.append(np.log10(np.clip(values, 1e-10, None)))  # Avoid log(0)

        # Normalize the distributions to match ASCII space
        x_min = min([np.min(dist) for dist in distributions])
        x_max = max([np.max(dist) for dist in distributions])
        x_values = np.linspace(x_min - (0.1 * abs(x_min)) , x_max + (0.1 * abs(x_max)), 10)

        # Generate ASCII plot using gnuplotlib
        plots = []
        for dist, label in zip(distributions, labels):
            y_values, _ = np.histogram(dist, bins=10, range=(x_min - (0.5 * abs(x_min)) , x_max + (0.5 * abs(x_max))), density=True)
            plots.append(
                (x_values, y_values, dict(legend=label, _with="lines"))
            )

        gp.plot(
            *plots,
            terminal="dumb 100,30",  # ASCII art with appropriate dimensions
            title="Chi-square Distributions (ASCII Art)",
            xlabel="Log(Chiscore)",
            ylabel="Density",
            set=["grid"],  # Apply a uniform grid
        )

    def plot_distributions(self):
        labels = []
        distributions = []
        for label, values in self.data.items():
            labels.append(label)
            distributions.append(values)

        plt.figure(figsize=(10, 6))
        colors = ['magenta', 'blue', 'green', 'orange']
        for idx, (dist, label) in enumerate(zip(distributions, labels)):
            sns.kdeplot(
                np.array(dist),
                fill=True,
                label=label,
                color=colors[idx % len(colors)],
                alpha=0.5
            )

        plt.xscale('log')
        plt.xlabel('Chiscore')
        plt.ylabel('Density')
        plt.title('Density Plot of Chi-Square Scores')
        plt.legend(title='Distributions')
        plt.grid(True, which="both", linestyle='--', linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        plt.savefig(self.png_file, dpi=300)
        plt.close()

if __name__ == '__main__':
    plotter = ChiSquarePlotter(json_file='tests/data/example_global_s1_pruned_score_dict.json', png_file='tests/data/example_global_s1_pruned_score.png')
    plotter.plot_distributions()
    plotter.plot_distributions_ascii()