from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd


class GeneExpressionData:
    """
    Parent class for handling gene expression data.
    """

    def __init__(self, file_path):
        """
        Read the expression data file into a pandas DataFrame.
        and set the first column as the index.
        """
        self.data = pd.read_csv(file_path, sep="\t")
        self.data.set_index(self.data.columns[0], inplace=True)

    def get_column_names(self):
        """
        Return the list of column names in the DataFrame.
        """
        return list(self.data.columns)


class RNASeqAnalysis(GeneExpressionData):
    """
    Child class for analyzing RNA-Seq gene expression data.
    """

    def compute_column_means(self):
        """
        Compute mean expression for each column in the DataFrame.
        """
        return self.data.mean()

    def plot_column_means(self, output_file):
        """
        Generate a bar plot of mean expression per column
        and save it as a figure to disk.
        """
        means = self.compute_column_means()
        plt.figure()
        means.plot(kind="bar")
        plt.title("Mean Gene Expression (Raw Data)")
        plt.ylabel("Mean expression")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()

    def filter_genes_by_threshold(self, threshold):
        """
        Filter genes based on a user-defined expression threshold.
        A gene is kept if its mean expression is above the threshold.
        """
        gene_means = self.data.mean(axis=1)
        filtered_data = self.data.loc[gene_means > threshold]
        return filtered_data

    def filter_genes_by_threshold_optimized(self, threshold):
        """
        Alternative implementation of the same mean-based gene filter.
        """
        return self.data[self.data.mean(axis=1) > threshold]

    def plot_filtered_column_means(self, threshold, output_file):
        """
        Plot and save mean expression per column
        for the filtered data.
        """
        filtered_data = self.filter_genes_by_threshold(threshold)
        filtered_means = filtered_data.mean()
        plt.figure()
        filtered_means.plot(kind="bar")
        plt.title(f"Mean Gene Expression (Filtered, Threshold={threshold})")
        plt.ylabel("Mean expression")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()


def main():
    """Generate plots from the coursework count matrix."""
    project_dir = Path(__file__).resolve().parent
    output_dir = project_dir / "figures"
    output_dir.mkdir(exist_ok=True)
    analysis = RNASeqAnalysis(project_dir / "data" / "GSE81082_count_matrix.txt")
    threshold = 1000
    filtered = analysis.filter_genes_by_threshold(threshold)
    analysis.plot_column_means(output_dir / "mean_expression_raw.png")
    analysis.plot_filtered_column_means(
        threshold, output_dir / "mean_expression_filtered.png"
    )
    print(f"Samples: {analysis.data.shape[1]}")
    print(f"Genes before filtering: {analysis.data.shape[0]}")
    print(f"Genes after filtering: {filtered.shape[0]}")
    print(f"Figures saved to: {output_dir}")


if __name__ == "__main__":
    main()
