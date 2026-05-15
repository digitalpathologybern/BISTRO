import upsetplot
import argparse
import matplotlib.pyplot as plt
import pandas as pd

def hvg_upset_plot(hvg_cat, top_percentage, output_plot, output_csv):

    if len(hvg_cat) == 0:
        raise ValueError('No csv files provided')

    for fraction in top_percentage:
        hvg_indicator = pd.DataFrame()
        methods = []
        for i, csv in enumerate(hvg_cat):
            normMethod = csv.split('/')[-1].split('.')[0].split('_')[-2]
            methods.append(normMethod)
            df = pd.read_csv(csv, index_col=0)
            # Make sure that the gene names do not contain - or / and substitute them by a .
            df.index = df.index.str.replace('-', '.')
            df.index = df.index.str.replace('/', '.')
            df.index = df.index.str.replace(' ', '.')
            df.index = df.index.str.replace(':', '.')
            df.index = df.index.str.replace('_', '.')

            ngenes = df.shape[0]

            # Sort the genes based on the dispersion values in descending order
            df = df.sort_values(by='dispersions', ascending=False)
            if i == 0:
                # Set index to the gene names
                hvg_indicator = pd.DataFrame(index=df.index)

                # Add a column normMethod to the dataframe and set to True the first 10% genes
                hvg_indicator[normMethod] = False
                hvg_indicator.loc[df.index[:int(ngenes*fraction)], normMethod] = True

            else:
                hvg_indicator[normMethod] = False
                
                # Set to True in hvg_indicator, the genes that are in the top 10% of the dispersion values in the df taking into account that the order may be different
                # for gene in df.index[:int(ngenes*0.1)]:
                hvg_indicator.loc[df.index[:int(ngenes*fraction)], normMethod] = True

        # Create a figure and plot the upset plot
        fig, ax = plt.subplots()
        # Remoev the axis from figure
        ax.axis('off')
        upsetplot.UpSet(upsetplot.from_indicators(methods, data=hvg_indicator), show_counts=True).plot(fig=fig)
        plt.savefig(output_plot  +'/top_'+str(fraction)+'_hvg_upset_plot.png')

        # Save the df as a csv
        hvg_indicator.to_csv(output_csv + '/top_'+str(fraction)+'_selected_hvg_per_method.csv')




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot upset plot for highly variable genes')
    parser.add_argument('--hvg_csv', nargs='+', default = [], help='List of paths to the csv with the genes and the dispersion values')
    parser.add_argument('--top_fraction', nargs='+', default = [0.1], type = float, help='Fraction of genes to be considered as highly variable')
    parser.add_argument('--output_plot', type = str, help='Output directory to save the plot')
    parser.add_argument('--output_csv', type = str, help='Output directory to save the csv')
    args = parser.parse_args()

    hvg_csv = args.hvg_csv
    top_fraction = args.top_fraction
    output_plot = args.output_plot
    output_csv = args.output_csv
    
    hvg_upset_plot(hvg_csv, top_fraction, output_plot, output_csv)