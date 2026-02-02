import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

class PhylogeneticBenchmark:
    """ Extract and compare results from CellPhy and IQ-TREE 3 for the cancer models. """

    def __init__(self, base_dir, genotype_model='GT10', output_dir=None):
        """
        Docstring for __init__
        
        :param self: Description
        :param base_dir: Base directory containing sim1, sim2 etc. folders.
        :param genotype_model: Genotype model used in the simulations - default is 'GT10'.
        """
        self.base_dir = Path(base_dir)
        self.genotype_model = genotype_model
        self.results = []

        # Set output directory
        if output_dir is None:
            self.output_dir = Path.cwd()
        else:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output files will be saved to: {self.output_dir}")

        # Define genotype states based on the model
        if genotype_model == 'GT10':
            self.genotypes = ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K']
        elif genotype_model == 'GT16':
            # TODO: check the order of genotypes in GT16. 
            self.genotypes = ['AA', 'CC', 'GG', 'TT', 
                            'AC', 'CA', 'AG', 'GA', 'AT', 'TA', 
                            'CG', 'GC', 'CT', 'TC', 'GT', 'TG']
        else:
            raise ValueError("Unsupported genotype model. Use 'GT10' or 'GT16'.")

    def extract_cellphy_likelihood(self, cellphy_log_file):
        """ Extract the log-likelihood from CellPhy log file. """
        try:
            with open(cellphy_log_file, 'r') as file:
                content = file.read()
                # Search for the final log-likelihood value
                match = re.search(r'Final LogLikelihood:\s+([-\d.]+)', content)
                if not match:
                    match = re.search(r'Final Log Likelihood:\s+([-\d.]+)', content)
                if match:
                    return float(match.group(1))
        except Exception as e:
            print(f"Error reading {cellphy_log_file}: {e}")
        return None
    
    def extract_cellphy_frequencies(self, log_file):
        """ Extract nucleotide frequencies from CellPhy model file. """
        try:
            with open(log_file, 'r') as file:
                content = file.read()
                # look for frequencies line
                freqs = {}

                pattern = r'Base frequencies.*?:\s+((?:[\d.eE+-]+\s*)+)'
                match = re.search(pattern, content)

                if match:
                    freq_str = match.group(1).strip()
                    freq_values = [float(x) for x in freq_str.split()]

                    if len(freq_values) == len(self.genotypes):
                        freqs = {self.genotypes[i]: freq_values[i] for i in range(len(self.genotypes))}
                        return freqs
        except Exception as e:
            print(f"Error reading {log_file}: {e}")
        return None
    
    def extract_cellphy_runtime(self, log_file):
        """ Extract runtime from CellPhz log file. """
        try:
            with open(log_file, 'r') as file:
                content = file.read()
                # Common patterns for runtime in RAxML-NG/ CellPhy logs
                # "Elapsed time: X seconds"
                pattern = r'Elapsed time:\s+([\d.]+)\s+seconds'
                match = re.search(pattern, content)
                if match:
                    return float(match.group(1))
        except Exception as e:
            print(f"Error reading {log_file}: {e}")
        return None
        
    def extract_iqtree_likelihood(self, iqtree_log_file):
        """ Extract the log-likelihood from IQ-TREE log file. """
        try:
            with open(iqtree_log_file, 'r') as file:
                content = file.read()
                # Search for the final log-likelihood value
                match = re.search(r'BEST SCORE FOUND :\s+([-\d.]+)', content)
                if not match:
                    match = re.search(r'Optimal log-likelihood:\s+([-\d.]+)', content)
                if match:
                    return float(match.group(1))
        except Exception as e:
            print(f"Error reading {iqtree_log_file}: {e}")
        return None
    
    def extract_iqtree_frequencies(self, iqtree_file):
        """Extract genotype frequencies from IQ-TREE .iqtree file"""
        try:
            with open(iqtree_file, 'r') as f:
                content = f.read()
                
            # IQ-TREE uses IUPAC codes in format (multi-line with indentation):
            # State frequencies: (estimated with maximum likelihood)
            #   pi(A) = 0.2893
            #   pi(C) = 0.1991
            #   pi(G) = 0.2006
            #   ...
            
            # Extract all pi(X) = value pairs (can span multiple lines)
            pattern = r'pi\(([A-Z])\)\s*=\s*([\d.eE+-]+)'
            matches = re.findall(pattern, content)
            
            if matches:
                freqs = {}
                for genotype, freq_str in matches:
                    # Convert IUPAC code to genotype notation
                    if genotype in self.genotypes:
                        freqs[genotype] = float(freq_str)
                
                if len(freqs) == len(self.genotypes):
                    return freqs
                else:
                    print(f"  Warning: Found {len(freqs)} frequencies but expected {len(self.genotypes)}")
                    print(f"  Found: {list(freqs.keys())}")
                        
        except Exception as e:
            print(f"Error reading {iqtree_file}: {e}")
        return None
    
    def extract_iqtree_runtime(self, iqtree_log_file):
        """ Extract runtime from IQ-TREE log file. """
        try:
            with open(iqtree_log_file, 'r') as file:
                content = file.read()
                # Common patterns for runtime in IQ-TREE logs
                # "Total wall-clock time used: X seconds"
                pattern = r'Total wall-clock time used:\s+([\d.]+)\s*(?:seconds?|sec)'
                match = re.search(pattern, content)
                if match:
                    return float(match.group(1))
        except Exception as e:
            print(f"Error reading {iqtree_log_file}: {e}")
        return None
    
    def calculate_frequency_distance(self, freq1, freq2):
        """ Calculate the Euclidean distance between two frequency distributions. """
        if freq1 is None or freq2 is None:
            return None
        
        vec1 = np.array([freq1.get(g, 0) for g in self.genotypes])
        vec2 = np.array([freq2.get(g, 0) for g in self.genotypes])

        return np.linalg.norm(vec1 - vec2)
    
    def process_datasets(self, sim_num, subsample):
        """
        Process datasets for a given simulation number and subsample. (e.g.  sim1/sim1.D0.00G0.00j250)

        :param sim_num: simulation numer, e.g. 'sim1'
        :param subsample: subsample identifier, e.g. 'D0.00G0.00j250'
        """
        sim_dir = self.base_dir / f"{sim_num}/{sim_num}.{subsample}"
        
        print(f"Looking for: {sim_dir}")

        if not sim_dir.exists():
            print(f"Directory {sim_dir} does not exist.")
            return
        cellphy_dir = sim_dir / "cellphy0.9.3"
        iqtree_dir = sim_dir / "iqtree3"

        print(f"CellPhy dir exists: {cellphy_dir.exists()}")
        print(f"IQ-TREE dir exists: {iqtree_dir.exists()}")

        if not cellphy_dir.exists() or not iqtree_dir.exists():
            print(f"Missing subdirectories of CellPhy or IQ-TREE in {sim_dir}.")
            return
        
        # Find all replicates
        cellphy_logs = sorted(cellphy_dir.glob("true_hap.*.raxml.log"))
        print(f"Found {len(cellphy_logs)} CellPhy log files.")

        if len(cellphy_logs) == 0:
            print(f"No CellPhy log files found in {cellphy_dir}.")
            print(f"Files in directory: {list(cellphy_dir.glob('*'))[:5]}")
            return
        
        for log_file in cellphy_logs:
            # Extract replicate number
            match = re.search(r'true_hap\.(\d+)\.', log_file.name)
            if not match:
                continue
            replicate = int(match.group(1))

            # CellPhy files
            cellphy_log = log_file
            cellphy_model = cellphy_dir / f"true_hap.{replicate:04d}.cellphy93.gt10.raxml.bestModel"

            # IQ-TREE files
            iqtree_log = iqtree_dir / f"true_hap.{replicate:04d}.iqtree3.gt10FO.log"
            iqtree_model = iqtree_dir / f"true_hap.{replicate:04d}.iqtree3.gt10FO.iqtree"

            # Extract data
            cellphy_logL = self.extract_cellphy_likelihood(cellphy_log)
            cellphy_freqs = self.extract_cellphy_frequencies(cellphy_log)
            #print(f'CellPhy frequencies: {cellphy_freqs}')
            cellphy_time = self.extract_cellphy_runtime(cellphy_log)

            iqtree_logL = self.extract_iqtree_likelihood(iqtree_log)
            iqtree_freqs = self.extract_iqtree_frequencies(iqtree_model)
            #print(f'IQ-TREE frequencies: {iqtree_freqs}')
            iqtree_time = self.extract_iqtree_runtime(iqtree_log)

            # Calculate metrics
            delta_logL = None
            if cellphy_logL is not None and iqtree_logL is not None:
                delta_logL = cellphy_logL - iqtree_logL
            
            freq_distance = self.calculate_frequency_distance(cellphy_freqs, iqtree_freqs)

            # Store results
            result = {
                'sim_num': sim_num,
                'subsample': subsample,
                'replicate': replicate,
                'cellphy_logL': cellphy_logL,
                'iqtree_logL': iqtree_logL,
                'delta_logL': delta_logL,
                'freq_distance': freq_distance,
                'cellphy_time': cellphy_time,
                'iqtree_time': iqtree_time
            }

            # Add individual frequencies
            if cellphy_freqs:
                for g in self.genotypes:
                    result[f'cellphy_freq_{g}'] = cellphy_freqs.get(g, None)
            if iqtree_freqs:
                for g in self.genotypes:
                    result[f'iqtree_freq_{g}'] = iqtree_freqs.get(g, None)

            self.results.append(result)

    def process_all(self, sim_num, rep_range=range(1, 101), subsamples=['D0.00G0.00j250', 'D0.00G0.00j500', 'D0.00G0.00j1000']):
        """
        Process all datasets

        :param sim_num: Simulation numbers to process (default: ['sim1', 'sim2'])
        :param rep_range: Replicate numbers to process (default: 1 to 100)
        :param subsamples: Subsample sizes to process
        """
        total = len(rep_range) * len(subsamples)
        count = 0

        for rep_num in rep_range:
            for subsample in subsamples:
                count += 1
                print(f"Processing {rep_num} {subsample} ({count}/{total})")
                self.process_datasets(sim_num,subsample)
        print(f"\nTotal results extracted> {len(self.results)}")

    def get_dataframes(self):
        """ Convert results to pandas DataFrame. """
        return pd.DataFrame(self.results)
        
    def save_results(self, filename="benchmark_results.csv"):
        """ Save results to a CSV file. """
        df = self.get_dataframes()
        output_path = self.output_dir / filename 
        print(f"\nSaving results to: {output_path.absolute()}")
        df.to_csv(output_path, index=False)

    def generate_summary_stats(self):
        """ Generate summary statistics for the results. """
        df = self.get_dataframes()

        if len(df) == 0:
            print("\n" + "="*60)
            print("WARNING: No data was extracted!")
            print("="*60)
            print("\nPlease check:")
            print("1. Base directory path is correct")
            print("2. Directory structure matches expected format")
            print("3. File naming conventions are correct")
            print("\nRun with a single dataset first to debug:")
            print("  benchmark.process_dataset(1, 'j250')")
            return
        
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60 + "\n")

        print(f"\nTotal datasets processed: {len(df)}")
        print(f"Datasets with valid likelihoods: {df['delta_logL'].notna().sum()}")
        print(f"Datasets with valid frequencies: {df['freq_distance'].notna().sum()}")

        print("\n--- Log-Likelihood Differences (CellPhy - IQ-TREE) ---")
        print(df['delta_logL'].describe())

        print("\n--- Frequency Euclidean Distances ---")
        print(df['freq_distance'].describe())

        print("\n--- By Subsample ---")
        for subsample in df['subsample'].unique():
            subset = df[df['subsample'] == subsample]
            print(f"\n{subsample}")
            print(f" Mean ΔlogL: {subset['delta_logL'].mean():.4f}")
            print(f" Std ΔlogL: {subset['delta_logL'].std():.4f}")
            print(f"Mean freq distance: {subset['freq_distance'].mean():.6f}")

    def plot_results(self, filename_prefix="benchmark_plot"):
        """Generate visualization plots"""
        df = self.get_dataframes()
        
        if len(df) == 0:
            print("\n✗ Cannot generate plots: No data available!")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Plot 1: Distribution of log-likelihood differences
        ax = axes[0, 0]
        df['delta_logL'].hist(bins=50, ax=ax, edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel('ΔlogL (CellPhy - IQ-TREE)')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Log-Likelihood Differences')
        
        # Plot 2: Log-likelihood differences by subsample size
        ax = axes[0, 1]
        df.boxplot(column='delta_logL', by='subsample', ax=ax)
        ax.set_xlabel('Subsample Size')
        ax.set_ylabel('ΔlogL (CellPhy - IQ-TREE)')
        ax.set_title('Log-Likelihood Differences by Subsample Size')
        plt.sca(ax)
        plt.xticks(rotation=0)
        
        # Plot 3: Frequency distance distribution
        ax = axes[1, 0]
        df['freq_distance'].hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel('Euclidean Distance between Frequencies')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Frequency Differences')
        
        # Plot 4: Frequency distance by subsample size
        ax = axes[1, 1]
        df.boxplot(column='freq_distance', by='subsample', ax=ax)
        ax.set_xlabel('Subsample Size')
        ax.set_ylabel('Frequency Distance')
        ax.set_title('Frequency Differences by Subsample Size')
        plt.sca(ax)
        plt.xticks(rotation=0)
        
        plt.tight_layout()
        output_path = self.output_dir / f'{filename_prefix}.png'
        print(f"\nSaving plots to: {output_path.absolute()}")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        # Verify the file was created
        if output_path.exists():
            file_size = output_path.stat().st_size
            print(f"✓ Plots saved successfully ({file_size:,} bytes)")
        else:
            print(f"✗ ERROR: Plot file was not created!")
        
        plt.close()
        return output_path

    def plot_likelihood_comparison(self, filename_prefix="likelihood_comparison.png"):
        """ Generate line plots comparing log-likelihoods between methods (CellPhy vs IQ-TREE). """
        df = self.get_dataframes()
        
        if len(df) == 0:
            print("\n✗ Cannot generate likelihood comparison plot: No data available!")
            return
        
        # Ensure output_dir is set
        if not hasattr(self, 'output_dir'):
            self.output_dir = Path.cwd()
            print(f"Warning: output_dir not set. Using current working directory: {self.output_dir}")

        # Get unique subsamples
        subsamples = sorted(df['subsample'].unique())
        n_subsamples = len(subsamples)

        # Create figure with subplots for each subsample
        fig, axes = plt.subplots(n_subsamples, 1, figsize=(14, 5 * n_subsamples))

        # If only one subsample, axes won't be listed
        if n_subsamples == 1:
            axes = [axes]

        for idx, subsample in enumerate(subsamples):
            ax = axes[idx]
            subset = df[df['subsample'] == subsample].sort_values('replicate')

            # Plot both methods
            ax.plot(subset['replicate'], subset['cellphy_logL'], label='CellPhy', marker='o', alpha=0.7, linewidth=2)
            ax.plot(subset['replicate'], subset['iqtree_logL'], label='IQ-TREE', marker='s', alpha=0.7, linewidth=2)
            ax.set_xlabel('Replicate Number')
            ax.set_ylabel('Log-Likelihood')
            ax.set_title(f'Log-Likelihood Comparison for Subsample: {subsample}')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_path = self.output_dir / filename_prefix
        print(f"\nSaving likelihood comparison plot to: {output_path.absolute()}")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
       
        # Verify the file was created
        if output_path.exists():
            file_size = output_path.stat().st_size
            print(f"✓ Likelihood comparison plot saved successfully ({file_size:,} bytes)")
        else:
            print(f"✗ ERROR: Likelihood comparison plot file was not created!")
        
        plt.close()
        return output_path

    def plot_runtime_comparison(self, filename_prefix="runtime_comparison.png"):
        """ Generate bar plots comparing runtimes between methods (CellPhy vs IQ-TREE). """
        df = self.get_dataframes()
        
        if len(df) == 0:
            print("\n✗ Cannot generate runtime comparison plot: No data available!")
            return
        
        # Ensure output_dir is set
        if not hasattr(self, 'output_dir'):
            self.output_dir = Path.cwd()
            print(f"Warning: output_dir not set. Using current working directory: {self.output_dir}")

        # Get unique subsamples
        subsamples = sorted(df['subsample'].unique())
        n_subsamples = len(subsamples)

        # Create figure with subplots for each subsample
        #fig, axes = plt.subplots(n_subsamples, 1, figsize=(14, 5 * n_subsamples))
        fig, axes = plt.subplots(2, n_subsamples, figsize=(7 * n_subsamples, 10))

        # If only one subsample, axes won't be listed
        if n_subsamples == 1:
            axes = axes.reshape(-1,1)

        for idx, subsample in enumerate(subsamples):
            subset = df[df['subsample'] == subsample].sort_values('replicate')

            # Top row: Line Plot of runtime across replicates
            ax = axes[0, idx]
            ax.plot(subset['replicate'], subset['cellphy_time'], label='CellPhy', marker='o', alpha=0.7, linewidth=2)
            ax.plot(subset['replicate'], subset['iqtree_time'], label='IQ-TREE', marker='s', alpha=0.7, linewidth=2)
            ax.set_xlabel('Replicate Number')
            ax.set_ylabel('Runtime (seconds)')
            ax.set_title(f'Runtime Comparison for Subsample: {subsample}')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # Bottom row: Box Plot of average runtime
            ax = axes[1, idx]
            runtime_data = [subset['cellphy_time'].dropna(), subset['iqtree_time'].dropna()]
            bp = ax.boxplot(runtime_data, labels=['CellPhy', 'IQ-TREE'], patch_artist=True)
            bp['boxes'][0].set_facecolor('lightblue')
            bp['boxes'][1].set_facecolor('lightgreen')
            ax.set_ylabel('Runtime (seconds)')
            ax.set_title(f'Runtime Distribution for Subsample: {subsample}')
            ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        output_path = self.output_dir / filename_prefix
        print(f"\nSaving runtime comparison plot to: {output_path.absolute()}")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

        # Verify the file was created
        if output_path.exists():
            file_size = output_path.stat().st_size
            print(f"✓ Runtime comparison plot saved successfully ({file_size:,} bytes)")
        else:
            print(f"✗ ERROR: Runtime comparison plot file was not created!")

        plt.close()
        return output_path

if __name__ == "__main__":
    # set the base directory
    base_dir = "/Users/u7875558/Documents/promotion/projects/cancer_models/benchmark_gt10/simulation_data_CellPhy/"
    # choose genotype model
    genotype_model = 'GT10'
    # set up the output directory
    output_dir = os.path.join(base_dir, "benchmark_results")
    os.makedirs(output_dir, exist_ok=True)
    #output_dir = base_dir / "benchmark_results"
    # initialise benchmark
    benchmark = PhylogeneticBenchmark(base_dir, genotype_model=genotype_model, output_dir=output_dir)

    # Process all datasets (or specify a subset for testing)
    # For testing, try benchmark.process_datasets(1, 'D0.00G0.00j250')
    #benchmark.process_all(sim_range=range(1, 11), subsamples=['D0.00G0.00j250', 'D0.00G0.00j500', 'D0.00G0.00j1000'])
    benchmark.process_datasets('sim2', 'D0.00G0.00')
    #benchmark.process_all(sim_num='sim1', rep_range=range(1, 101), subsamples=['D0.00G0.00j250', 'D0.00G0.00j500', 'D0.00G0.00j1000'])
    # save results
    benchmark.save_results("sim2_benchmark_results_gt10FO.csv")
    
    # generate summary statistics
    benchmark.generate_summary_stats()
    
    # generate plots
    benchmark.plot_results(filename_prefix='sim2_benchmark_gt10FO_plots')

    # generate likelihood comparison plots
    benchmark.plot_likelihood_comparison(filename_prefix='sim2_likelihood_comparison_gt10FO.png')

    # generate runtime comparison plots
    benchmark.plot_runtime_comparison(filename_prefix='sim2_runtime_comparison_gt10FO.png')

    # Access the DataFrame 
    df = benchmark.get_dataframes()
    print("\nFirst few rows:")
    print(df.head())