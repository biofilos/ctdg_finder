import os
import sys
from glob import glob

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# from matplotlib import patches

"""
1 argument: ctdg database directory
2 argument: directory with merged ctdg clusters
3 argument: species name
"""
# Update global matplotlib parameters
params = {"axes.grid": False}
mpl.rcParams.update(params)

class CgpOut:
    def __init__(self, in_table, genomes, all_genes_tab, sp):
        def get_chrom(sp_cl):
            """
            Extract chromosome name from a tuple containing species and cluster
            :param sp_cl: list = [sp, cluster]
            :return: chromosome name
            """
            return self.only_clustered.loc[(self.only_clustered['species'] == sp_cl[0]) &
                                           (self.only_clustered['cluster'] == sp_cl[1]),
                                           "chromosome"].values[0]

        def get_chrom_len(sp_chrom):
            """
            Get chromosome length
            :param sp_chrom: list([species, chromosome])
            :return: chromosome length
            """
            return self.genomes.loc[(self.genomes['species'] == sp_chrom[0]) &
                                    (self.genomes['chromosome'] == sp_chrom[1]), "length"].values[0]

        def get_num_genes(chrom_cluster):
            """
            Get the number of genes for a specific category (clustered, single-copy, etc)
            :param chrom_cluster: pd.Series with two fields: chromosome, status
            :return: number of genes in that category for a cluster
            """
            chrom, status = chrom_cluster.values
            try:
                return genes_in_cluster.loc[(chrom, status), 'acc']
            except KeyError:
                return 0

        self.sp_genes = all_genes_tab.loc[all_genes_tab['species'] == sp]
        self.cgps = in_table.loc[in_table['species'] == sp]
        self.genomes = genomes
        # annotate clustered genes and non-clustered genes
        self.cgps.loc[:, 'order'] = self.cgps['order'].astype(str)
        # Annotate types of genes in a cluster
        # Single copy genes (in chromosome)
        self.cgps.loc[self.cgps['order'] == "0.0", "status"] = "single_copy"
        # Genes rejected by MeanShift
        self.cgps.loc[self.cgps['order'] == "na_ms", "status"] = "rejected_ms"
        # Genes rejected by the statistical assessment
        self.cgps.loc[self.cgps['order'] == "na_95", "status"] = "rejected_95"
        # Clustered genes
        self.cgps.loc[self.cgps['status'].isnull(), "status"] = "clustered"
        # Extract genes that are clustered
        self.htop
        only_clustered = self.cgps.loc[self.cgps['status'] == "clustered"]
        # Calculate number of genes in cluster from each category
        genes_in_cluster = self.cgps.groupby(['chromosome', 'cluster', 'status']).count().reset_index()
        genes_in_cluster.set_index(['cluster', 'status'], inplace=True)
        first_in_cluster = self.only_clustered.groupby(["species", "cluster"]).min()["start"]
        last_in_cluster = self.only_clustered.groupby(["species", "cluster"]).max()["end"]
        self.clusters = pd.concat([first_in_cluster, last_in_cluster], 1)  # type: pd.DataFrame
        # Annotate chromosome names
        self.clusters.loc[:, "chromosome"] = self.clusters.index.map(get_chrom)
        self.clusters.loc[:, "chromosome"] = self.clusters['chromosome'].astype(str)
        self.clusters.reset_index(inplace=True)
        # Get chromosome lengths
        chromosome_lengths = self.clusters[["species", "chromosome"]].apply(get_chrom_len, 1)
        self.clusters.loc[:, "chromosome_length"] = chromosome_lengths
        self.chromosomes = set(self.clusters['chromosome'])
        # Add number of genes for each category in each cluster
        for status in set(self.cgps['status']):
            # Insert dummy column so that I can pass both arguments as one
            self.clusters.loc[:, 'dummy'] = status
            self.clusters.loc[:, status] = self.clusters[['cluster', 'dummy']].apply(get_num_genes, 1)
        del self.clusters['dummy']
        # Calculate length of each cluster
        self.clusters.loc[:, "cluster_length"] = self.clusters['end'] - self.clusters['start']
        lengths_100 = self.clusters['cluster_length'] * 100
        self.clusters.loc[:, "perc_region_clustered"] = lengths_100 / self.clusters['chromosome_length']
        sp_genes_dict = self.sp_genes.groupby('chromosome').count()['acc'].to_dict()
        self.clusters.loc[:, "total_genes"] = self.clusters['chromosome'].map(lambda x: sp_genes_dict[x])
        genes_100 = self.clusters['clustered'] * 100
        self.clusters.loc[:, "perc_clustered_genes"] = genes_100 / self.clusters['total_genes']

    def plot_chromosome_clusters(self, chrom):
        chrom_table = self.clusters.loc[self.clusters['chromosome'] == chrom]
        gene_numbers = self.only_clustered.groupby(["species", "cluster"]).count()['start'].to_dict()
        # Initialize figure
        fig = plt.figure()
        # Set axis for plotting the chromosome information
        # The argument is a list (x_position,y_position, width, height) as proportions of the image
        ax = fig.add_axes([0, 0, 1, 0.15])
        plt.title("Chromosome {}".format(chrom))
        num_genes_list = []
        rectangles = []
        # Plot each cluster in a chromosome
        for cluster in chrom_table.index:
            cl_tab = chrom_table.loc[cluster]
            cl_start = cl_tab["start"]
            cl_end = cl_tab["end"]
            num_genes = gene_numbers[(cl_tab["species"], cl_tab["cluster"])]
            # Append the number of genes in the cluster to a list
            num_genes_list.append(num_genes)
            # Add a rectangle for each cluster
            # Set the height of the bar as the number of genes
            #     rectangles.append(patches.Rectangle((cl_start,0), cl_end-cl_start,num_genes))
            rectangles.append(mpl.patches.Rectangle((cl_start, 0), cl_end - cl_start, 1))
        # Calculate maximum and minimum number of genes (for the color bar)
        max_genes = max(num_genes_list)
        min_genes = min(num_genes_list)
        # Set the color palette
        #cmap = sns.color_palette("Blues_d", max_genes)
        cmap = sns.cubehelix_palette(max_genes, start=0, rot=-0.5,gamma=0.7,hue=0.5,light=0.7,dark=0)
        cmap_legend = mpl.colors.ListedColormap(cmap)
        # Pick colors from the palette for each cluster
        colors = [cmap[x - 1] for x in num_genes_list]

        p = mpl.collections.PatchCollection(rectangles, match_original=False,
                                            facecolors=colors, edgecolor='none')
        ax.add_collection(p)
        # With the height of the bar as the number of genes
        # plt.ylim(0,max_genes)
        plt.ylim(0, 1)
        # plt.yticks(None)
        plt.setp(ax.get_yticklabels(), visible=False)

        beginning = chrom_table.min()['start']
        finale = chrom_table.max()['end']
        plt.xlim(beginning, finale)
        ax_color = fig.add_axes([0, 0.3, 1, .1])
        norm = mpl.colors.Normalize(vmin=min_genes, vmax=max_genes)
        color_bar = mpl.colorbar.ColorbarBase(ax_color,
                                                     cmap=cmap_legend,
                                                     norm=norm,
                                                     orientation="horizontal")
        color_bar.set_label("Genes per CTDG")
        return fig

    def plot_all_chromosome_clusters(self, chrom_set):
        gene_numbers = self.only_clustered.groupby(["species", "cluster"]).count()['start'].to_dict()
        max_genes = max(gene_numbers.values())
        min_genes = 0
        # Initialize figure
        fig = plt.figure(figsize=(8, 8))
        # Initialize plot coordinates
        ax_x = 0
        ax_y = False
        num_chroms = len(chrom_set)
        # If there are odd number of chromosomes, add one
        if num_chroms % 2 != 0:
            num_chroms += 1
            skip_even = True
        else:
            skip_even = False
        rows = int(num_chroms / 2) + 1
        # create grid of plots
        gs = mpl.gridspec.GridSpec(rows, 2)
        # Plot each chromosome
        for chrom in chrom_set:
            # Filter chromosome table
            chrom_table = self.clusters.loc[self.clusters['chromosome'] == chrom]
            # Place plot coordinates
            ax = plt.subplot(gs[ax_x, int(ax_y)])

            # Set plot title
            #             ax.set_title("Chromosome {}".format(chrom))
            #             gs.update(hspace=.5)
            # Initialize number of genes and patches list
            num_genes_list = []
            rectangles = []
            # Plot each cluster in a chromosome
            for cluster in chrom_table.index:
                cl_tab = chrom_table.loc[cluster]
                cl_start = cl_tab["start"]
                cl_end = cl_tab["end"]
                num_genes = gene_numbers[(cl_tab["species"], cl_tab["cluster"])]
                # Append the number of genes in the cluster to a list
                num_genes_list.append(num_genes)
                # Add a rectangle for each cluster
                # Set the height of the bar as the number of genes
                #     rectangles.append(patches.Rectangle((cl_start,0), cl_end-cl_start,num_genes))
                rectangles.append(mpl.patches.Rectangle((cl_start, 0), cl_end - cl_start, 1))
            ax.text(chrom_table['chromosome_length'].values[0] / 2, 0.5,
                    "Chromosome {}".format(chrom),
                    va='center', ha='center', alpha=1, color='gray')
            ax.tick_params(axis='both', which='major', labelsize=5, pad=0)
            ax.tick_params(axis='both', which='minor', labelsize=5, pad=0)
            ax.xaxis.get_offset_text().set_fontsize(5)
            #             ax.xaxis.get_offset_text().set_visible(False)

            # Calculate maximum and minimum number of genes (for the color bar)
            # Set the color palette
            cmap = sns.cubehelix_palette(max_genes, start=0, rot=-0.5,gamma=0.7,hue=0.5,light=0.7,dark=0)
            #cmap = sns.color_palette("Blues_d", max_genes)
            cmap_legend = mpl.colors.ListedColormap(cmap)
            # Pick colors from the palette for each cluster
            colors = [cmap[x - 1] for x in num_genes_list]

            p = mpl.collections.PatchCollection(rectangles, match_original=False,
                                                facecolors=colors, edgecolor='none')
            ax.add_collection(p)
            # With the height of the bar as the number of genes
            # plt.ylim(0,max_genes)
            ax.set_ylim(0, 1)
            plt.setp(ax.get_yticklabels(), visible=False)
            # Plot from the beginning to the end of each chromosome
            # Plot for the first to the last clusters is commented out
            beginning = 0  # chrom_table.min()['start']
            finale = chrom_table['chromosome_length'].values[0]  # chrom_table.max()['end']
            ax.set_xlim(beginning, finale)
            # If plot y-coordinate is 1, add one to the plot x-coordinate
            #             print(ax_x, rows)

            ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))

            if ax_x >= rows - 2 or (ax_x == rows - 3 and skip_even and ax_y):
                ax.set_xlabel("Genomic coordinates", fontsize=5,
                              horizontalalignment='center')
            # ax.xaxis.set_tick_params(pad=-1)
            if ax_y:
                ax_x += 1
            # Switch the plot y-coordinate between 1 and 0 
            ax_y = not ax_y

        if skip_even:
            ax_x += 1
        ax_color = plt.subplot(gs[ax_x, :])
        norm = mpl.colors.Normalize(vmin=min_genes, vmax=max_genes)
        color_bar = mpl.colorbar.ColorbarBase(ax_color,
                                              cmap=cmap_legend,
                                              norm=norm,
                                              orientation="horizontal")
        color_bar.set_label("Genes per CTDG")
        #         fig.tight_layout()
        gs.update(hspace=1)
        return fig
