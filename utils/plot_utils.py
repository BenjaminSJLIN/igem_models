import matplotlib.pyplot as plt
import seaborn as sns
import os

class Plotter:
    def __init__(self, output_dir):
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
    def plot_time_series(self, data, y_col, title, filename):
        """
        Plot time-series dynamics for different Arabinose concentrations.
        """
        plt.figure(figsize=(10, 6))
        sns.lineplot(
            data=data, 
            x='Time', 
            y=y_col, 
            hue='Ara_Concentration', 
            palette='viridis'
        )
        plt.title(title)
        plt.xlabel('Time (minutes)')
        plt.ylabel(y_col + ' Concentration')
        plt.grid(True, alpha=0.3)
        
        # Make the legend nicer
        plt.legend(title='[Ara]', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filepath}")

    def plot_dose_response(self, data, y_col, title, filename):
        """
        Plot steady-state dose-response curve.
        Assumes the last time point in the data is the steady state.
        """
        # Get the final time point for each concentration
        steady_state_data = data.groupby('Ara_Concentration').last().reset_index()
        
        plt.figure(figsize=(8, 6))
        plt.plot(
            steady_state_data['Ara_Concentration'], 
            steady_state_data[y_col], 
            marker='o', 
            linestyle='-', 
            color='royalblue',
            linewidth=2,
            markersize=8
        )
        plt.title(title)
        plt.xlabel('Arabinose Concentration [Ara]')
        plt.ylabel(f'Steady State {y_col} Concentration')
        plt.grid(True, alpha=0.3)
        
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filepath}")
