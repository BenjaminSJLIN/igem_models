import matplotlib.pyplot as plt
import seaborn as sns
import os

class Plotter:
    def __init__(self, output_dir):
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
    def plot_line(self, data, x_col, y_col, hue_col, title, filename, legend_title=None, x_label=None, log_x=False):
        """
        Plot a line graph grouped by a specific column (hue_col).
        """
        plt.figure(figsize=(10, 6))
        sns.lineplot(
            data=data, 
            x=x_col, 
            y=y_col, 
            hue=hue_col, 
            palette='viridis'
        )
        plt.title(title)
        
        if x_label is None:
            x_label = x_col
        plt.xlabel(x_label)
        plt.ylabel(y_col)
        
        if log_x:
            plt.xscale('log')
            
        plt.grid(True, alpha=0.3)
        
        # Make the legend nicer
        if legend_title is None:
            legend_title = hue_col
        plt.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filepath}")

    def plot_dose_response(self, data, x_col, y_col, title, filename, x_label=None, log_x=False):
        """
        Plot steady-state dose-response curve.
        Assumes the last time point in the data is the steady state, OR
        if 'Time' doesn't exist, just plots the data directly.
        """
        # Get the final time point for each concentration if 'Time' is present
        if 'Time' in data.columns:
            plot_data = data.groupby(x_col).last().reset_index()
        else:
            plot_data = data
            
        plt.figure(figsize=(8, 6))
        plt.plot(
            plot_data[x_col], 
            plot_data[y_col], 
            marker='o', 
            linestyle='-', 
            color='royalblue',
            linewidth=2,
            markersize=8
        )
        plt.title(title)
        
        if x_label is None:
            x_label = x_col
        plt.xlabel(x_label)
        plt.ylabel(f'Steady State {y_col}')
        
        if log_x:
            plt.xscale('log')
            
        plt.grid(True, alpha=0.3)
        
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filepath}")
