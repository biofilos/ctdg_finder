import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

in_genomes = sys.argv[1]
in_genomes_tab = pd.read_csv(in_genomes)

in_file = sys.argv[2]
in_tab = pd.read_csv(in_file)


class CgpOut:
    def __init__(self, in_table, genomes):
        self.cgps = in_table
        self.genomes = genomes
        # annotate clustered genes and non-clustered genes
        self.cgps.loc[self.cgps['order'] == "0.0", "status"] = "single_copy"
        self.cgps.loc[self.cgps['order'] == "na_ms", "status"] = "rejected_ms"
        self.cgps.loc[self.cgps['order'] == "na_95", "status"] = "rejected_95"
        self.cgps.loc[self.cgps['status'].isnull(), "status"] = "clustered"
        self.only_clustered = self.cgps.loc[self.cgps['status'] == "clustered"]
        self.clusters = pd.concat([self.only_clustered.groupby(["species", "cluster"]).min()["start"],
                                   self.only_clustered.groupby(["species", "cluster"]).max()["end"]], 1)

        def get_chrom(sp_cl):
            """
            Extract chromosome name from a tuple containing species and cluster
            :param sp_cl: list = [sp, cluster]
            :return: chromosome name
            """
            return self.only_clustered.loc[(self.only_clustered['species'] == sp_cl[0]) &
                                           (self.only_clustered['cluster'] == sp_cl[1]), "chromosome"].values[0]

        def get_chrom_len(sp_chrom):
            return self.genomes.loc[(self.genomes['species'] == sp_chrom[0]) &
                                    (self.genomes['chromosome'] == sp_chrom[1]), "length"].values[0]

        self.clusters.loc[:, "chromosome"] = self.clusters.index.map(get_chrom)
        self.clusters.reset_index(inplace=True)

        chromosome_lengths = self.clusters[["species", "chromosome"]].apply(get_chrom_len, 1)
        self.clusters.loc[:, "chromosome_length"] = chromosome_lengths

    def plot_cluster(self, sp, chrom):
        """
        Plot all the clusters in a chromosome
        :param sp: species name
        :param chrom: chromosome
        :return: None (save svg figure)
        """

        chrom_table = self.clusters.loc[(self.clusters['species'] == sp) &
                                        (self.clusters['chromosome'] == chrom)]





cgp = CgpOut(in_tab, in_genomes_tab)
