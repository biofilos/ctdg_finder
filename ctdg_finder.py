import argparse
import json
import os
import shutil
import tarfile
from concurrent import futures
from functools import partial
from glob import glob

import numpy as np
import pandas as pd
from sklearn.cluster import MeanShift

pd.options.mode.chained_assignment = None


# Define class CTDG
class CtdgConfig:
    def __init__(self, out_dir, db, hmmer_samples, sp):
        self.sp = sp
        self.hmmer_samples = hmmer_samples
        self.db = db
        self.out_dir = out_dir.strip('/')
        self.genome_file = "{}/chromosomes.csv".format(db)
        self.all_genes_file = "{}/genes_parsed.csv".format(db)
        # self.all_genes_fasta = "{}/all_seqs.fa".format(db)
        self.genomes = None
        self.all_genes = None
        self.hmmer_rows = None
        self.ms_result = None
        self.family_numbers_list = None
        # self.chrom_wide = []
        self.genome_wide = []

    def check_arguments(self):
        """
        Perform assertion statements to check if all the files needed
        for the analysis are present
        :return: None
        """
        # Check that the directory with the database exists
        assert os.path.exists(self.db), "directory {} with database does not exist".format(self.db)
        #db_path = '{}/pfam/pfam.hmm'.format(self.db)
        db_path = '{}/panther/panther.hmm'.format(self.db)
        # Check that the blast database files exists
        hmmer_files = sum([os.path.exists(db_path + x) for x in ['', '.h3f', '.h3i', '.h3m', '.h3p']])
        assert hmmer_files == 5, "Some of the HMMer database files with prefix pfam do not exist"
        # Check that the chromosomes information table exists
        assert os.path.exists(self.genome_file), "Chromosome information file" + \
                                                 "{}/chromosomes.csv does not exist".format(self.genome_file)
        # Check that the annotation table exists
        assert os.path.exists(self.genome_file), "Annotation file" + \
                                                 "{}/genes_parsed.csv does not exist".format(self.genome_file)

        # Check that number of blast samples has been set
        assert self.hmmer_samples, "Number of Hmmer samples was not specified"
        # Make sure that the output directory is not the directory where the database is or the directory where
        # the CTDGFinder resides
        banned_dirs = [os.getcwd(), self.db, '.']
        assert self.out_dir not in banned_dirs, "output direcrory can't be the ctdgfinder directory" + \
                                                " or the directory where the ctdg database is stored"

    def init_tables(self):
        """
        Initialize genome and gene annotation tables
        :return: None
        """
        genomes = pd.read_csv(self.genome_file)
        genomes.dropna(inplace=True)
        genomes.loc[:, 'chromosome'] = genomes['chromosome'].astype(str)
        all_genes = pd.read_csv(self.all_genes_file)
        all_genes.loc[:, 'chromosome'] = all_genes['chromosome'].astype(str)
        all_genes.set_index('acc', inplace=True)
        # self.all_genes = self.all_genes.loc[~(self.all_genes['chromosome'].isnull()) &
        #                                     (self.all_genes['chromosome'] != "NA")]

        # Filter species subset if necessary
        if self.sp:
            all_genes = all_genes.loc[all_genes["species"].isin(self.sp)]
            genomes = genomes.loc[genomes["species"].isin(self.sp)]
        self.genomes = genomes
        self.all_genes = all_genes


class CtdgRun:
    def __init__(self, ctdg, args):
        # From config tables
        self.genes = ctdg.all_genes
        self.genomes = ctdg.genomes
        self.out = ctdg.out_dir
        self.samples = ctdg.hmmer_samples
        # From arguments
        self.cpu = args.cpu - 1
        self.ref_seq = args.ref_seq
        self.ref_pfam = args.ref_pfam
        self.name = args.name_family
        self.db = args.db
        self.selected = None

    def create_folders(self):
        """
        Create the directories where all the
        results will be saved
        """
        if not os.path.exists(self.name):
            os.makedirs(self.name)
            folders = ['intermediates', 'report']
            for sub_dir in folders:
                target = self.name + '/' + sub_dir
                if not os.path.exists(target):
                    os.makedirs(target)
        print("Directory tree for {} was created.".format(self.name))

    def remove_if_unfinished(self):
        """
        Remove unfinished analysis
        :return: None
        """
        # Do not run analysis if final output is present
        assert not os.path.exists("{}/{}".format(self.out, self.name)), "Results for {} exist.".format(self.name)

        if os.path.exists(self.name):
            print("Partial results for {} were found, removing them".format(self.name))
            shutil.rmtree(self.name)

    def hmmscan(self, evalue=1e-3):
        out_file = "{0}/intermediates/{0}.pre_hmmer".format(self.name)
        dom_out = out_file + "_dom"
        long_out = out_file + "_long"
        #pfam_hmm = self.db + "/pfam/pfam.hmm"
        pfam_hmm = self.db + "/panther/panther.hmm"
        if self.ref_seq and not self.ref_pfam:
            print("Running HMMscan")
            cmd = "hmmscan -o {} --tblout {} --domtblout {} -E {} --cpu {} {} {}".format(long_out, out_file,
                                                                                         dom_out,
                                                                                         evalue, self.cpu,
                                                                                         pfam_hmm, self.ref_seq)
            # Run HMMerScan
            os.system(cmd)

            # Parse the Pfam accessions found
            pfams = []
            for line in open(out_file):
                if not line.startswith("#"):
                    data = [x for x in line.split(" ") if x != ""]
                    pfams.append(data[1])
        else:
            pfams = [self.ref_pfam]
        # Load Pfam Dictionary
        pfam_dict = json.load(open("{}/hmmer/pfam.json".format(self.db)))
        # Filter gene table with the genes that contain the Pfam domains found by HMMscan
        pfam_genes = []
        for pfam in pfams:
            if pfam in pfam_dict:
                pfam_genes += pfam_dict[pfam]
        # Remove duplicates
        pfam_selected = set(pfam_genes)
        # Subset table
        selected = self.genes.loc[pfam_selected]
        # Save clusters in list, so that they can be processed in parallel
        selected_list = []
        sp_chroms = selected.groupby(["species", "chromosome"]).groups
        for sp, chrom in sp_chroms:
            try:
                bandwidth = self.genomes.loc[(self.genomes["species"] == sp) &
                                             (self.genomes["chromosome"] == chrom),
                                             "bandwidth"].values[0]
            except IndexError:
                print(sp, chrom)
                # return []
            selected.loc[sp_chroms[(sp, chrom)], "bandwidth"] = bandwidth
            selected.loc[sp_chroms[(sp, chrom)], "cluster"] = "{}_{}".format(self.name, chrom)
            selected_list.append(selected.loc[sp_chroms[(sp, chrom)]])
        return selected_list

    @property
    def numbers(self):
        groups = self.selected.groupby(["species", "chromosome", "cluster"])
        starts = groups.min()["start"]
        ends = groups.max()["end"]
        dups = groups.count()["start"]
        numbers_df = pd.concat([dups, starts, ends], 1)
        numbers_df.columns = ["gene_duplicates", "start", "end"]
        numbers_df = numbers_df.loc[numbers_df["gene_duplicates"] > 1]
        numbers_df["length"] = numbers_df["end"] - numbers_df["start"]
        numbers_df.reset_index(inplace=True)
        return numbers_df


def meanshift(sp_table):
    """
    Performs meanshift clustering on the start coordinates of a table
    Importantly, all the genes in the input are considered to be on the same chromosome
    :param sp_table: pd.DataFrame of gene coordinates and bandwidth parameter
    for that chromosome
    :return: annotated table
    """
    # Order genes by start coordinate
    sp_table.sort_values("start", inplace=True)
    bandwidth = sp_table["bandwidth"].values[0]
    if bandwidth == 0:
        bandwidth = 1
    # If there is only one gene, it is single copy in the chromosome
    if type(sp_table) == pd.Series or sp_table.shape[0] == 1:
        sp_table.loc[:, "cluster"] = 0
        sp_table.loc[:, "order"] = 0
    else:
        gene_starts = np.array([(x, y) for x, y in zip(sp_table.start.values.astype(int),
                                                       np.zeros(len(sp_table)))])
        mean_shift = MeanShift(bandwidth=bandwidth)
        try:
            mean_shift.fit(gene_starts)
        except ValueError:
            print("BW: {}".format(bandwidth))
        labels = mean_shift.labels_
        # Label clusters
        sp_table.loc[:, "cluster"] = sp_table["cluster"] + ["_{}".format(x) for x in labels]
        # Set the order of each gene in the cluster
        cluster_groups = sp_table.groupby("cluster").groups
        for cluster in cluster_groups:
            # Get indices of the genes in a cluster
            cluster_ix = cluster_groups[cluster]
            genes_per_cluster = len(cluster_ix)
            if genes_per_cluster == 1:
                elements = ["na_ms"]
            else:
                elements = list(range(1, genes_per_cluster + 1))
            sp_table.loc[cluster_ix, "order"] = elements
    del sp_table["bandwidth"]
    return sp_table


def parallel_meanshift(hmm_df, cpu):
    """
    Run the MeanShift step in parallel
    :param hmm_df: parsed hmmscan table
    :param cpu: number of CPUs to use
    :return: hmmscan table with annotated cluster candidates
    """
    with futures.ProcessPoolExecutor(cpu) as pool:
        ms = pool.map(meanshift, hmm_df)
    return pd.concat(ms)


def sample_region(sample_chrom, record, chromosomes, genes, db):
    """
    Returns a tuple with chromosome, start and end
    :param record: list of records. species should be in position 1, chromosome in position 2, length
    in position 7
    :param genomes: Genomes annotation dataframe
    :return: tuple
    """
    sp = record[1]
    length = record[7]
    # Select from chromosomes at least the same size of the cluster
    # Get sample list of chromosomes
    # sample_chrom = np.random.choice(chromosomes.chromosome.values)
    sample_chrom_length = chromosomes.loc[chromosomes["chromosome"] == sample_chrom, "length"]
    sample_start = np.random.randint(1, sample_chrom_length - length)
    sample_end = sample_start + length

    sample_genes = genes.loc[(genes["chromosome"] == sample_chrom) &
                             (genes["start"] < sample_end) &
                             (genes["end"] > sample_start)].index
    json_path = "{}/hmmers/{}_{}.json".format(db, sp, sample_chrom)
    chrom_dict = json.load(open(json_path))
    max_dups_list = [len(set(chrom_dict[x]).intersection(set(sample_genes))) for x in sample_genes]
    if not max_dups_list:
        max_dups = 0
    else:
        max_dups = max(max_dups_list)
    # If there is only one gene, there are no duplicates
    if max_dups == 1:
        max_dups = 0
    # print("Chromosome {} ({} - {}): {}".format(sample_chrom, sample_start, sample_end, max_dups))
    return max_dups


def sample_record(record, ctdg_obj):
    """
    Calculates the sample distribution of duplicates for one record and returns the percentile 95
    :param record: tuple, where position 1: species, position 3: cluster, position 7: length of cluster, and
    position 8: the path to the chromosome JSON file
    :param ctdg-obj: CTDG object containing genes annotation, genome annotation, and db path
    :return: dictionary of the form: {(species,cluster): percentile 95
    """
    genes = ctdg_obj.genes
    genomes = ctdg_obj.genomes
    samples = ctdg_obj.samples
    db = ctdg_obj.db
    sp = record[1]
    # chrom = record[2]
    cluster = record[3]
    length = record[7]
    genomes = genomes.loc[(genomes["species"] == sp) &
                          (genomes["length"] >= length)]
    # Get array of chromosomes to sample
    sample_chroms = np.random.choice(genomes.chromosome, samples)
    genes = genes.loc[genes["species"] == sp]
    fill_sample_fx = partial(sample_region, record=record,
                             chromosomes=genomes, genes=genes, db=db)
    # Try with ProcessPoolExecutor
    with futures.ThreadPoolExecutor(ctdg_obj.cpu) as pool:
        max_dups = pool.map(fill_sample_fx, sample_chroms)
    percentile_95 = np.round(np.percentile(list(max_dups), 95), 4)
    return sp, cluster, percentile_95


def sample_table(ctdg_obj):
    """
    Annotate the cluster candidate summary table with the percentile 95 column
    :param numbers: cluster candidate table
    :param ctdg_obj: CTDG object
    :return: annotated table
    """
    numbers = ctdg_obj.numbers
    num_clusters = numbers.shape[0]
    if not num_clusters:
        return None
    print("Analyzing {} cluster candidates".format(num_clusters))

    longest_sp = max([len(x) for x in numbers.species.unique()]) + 2
    longest_cluster = max([len(x) for x in numbers.cluster.unique()]) + 2
    print_msg = "{:<{sp}} {:<{cluster}} {:<10} {}{}"
    print(print_msg.format("Species", "Cluster", "Duplicates", "P95", "",
                           sp=longest_sp, cluster=longest_cluster))
    for record in numbers.to_records():
        sp, cluster, p95 = sample_record(record, ctdg_obj)
        cluster_duplicates = record[4]
        if cluster_duplicates >= p95:
            msg = "*"
        else:
            msg = ""
        print(print_msg.format(sp, cluster, cluster_duplicates, p95, msg,
                               sp=longest_sp, cluster=longest_cluster))
        numbers.loc[(numbers["species"] == sp) & (numbers["cluster"] == cluster), "p_95"] = p95
    return numbers


def clean_p95(numbers):
    """
    Annotate genes in cluster candidates with less duplicates than
    the percentile 95 sampling, and save data
    :return:
    """
    if numbers.shape[0] > 0:
        for_removal_df = numbers.loc[numbers["p_95"] > numbers["gene_duplicates"]]
        # Clusters that did not pass the percentile 95 threshold
        for_removal_ix = for_removal_df.set_index(["species", "chromosome", "cluster"]).index
        genes = analysis.selected.reset_index().set_index(["species", "chromosome", "cluster"])
        # Annotate genes in clusters under the percentile 95 threshold
        genes.loc[for_removal_ix, "order"] = "na_p95"
        genes.reset_index(inplace=True)
        genes.set_index("acc", inplace=True)
        # Order columns
        genes.columns = ["species", "chromosome", "symbol", "start", "end", "strand",
                         "length", "cluster", "order"]
        # Save genes data
        genes_out_name = "{0}/report/{0}_genes.csv".format(analysis.name)
        genes.to_csv(genes_out_name)
        # Save only data of genes in clusters
        genes_clean = genes.loc[~(genes["order"].isin(["na_p95", "na_ms", 0, "0"]))]
        genes_clean.to_csv(genes_out_name.replace("genes", "genes_clean"))
        # Save clusters summary information

        numbers.set_index(["species", "chromosome", "cluster"], inplace=True)
        # Order columns
        numbers.columns = ["gene_duplicates", "start", "end", "length", "p_95"]
        numbers_out = genes_out_name.replace("genes", "numbers")
        numbers.to_csv(numbers_out)

        # Save only clusters that passed the percentile 95 sampling
        numbers_clean = numbers.loc[~(numbers.index.isin(for_removal_ix))]
        numbers_clean.to_csv(numbers_out.replace("numbers", "numbers_clean"))


def save_results(ctdg_object):
    # Compress intermediates
    out = ctdg_object.out
    fam_name = ctdg_object.name
    tar_name = "{}/intermediates/intermediates.tgz".format(fam_name)
    tar_file = tarfile.open(tar_name, "x:gz")
    for intermediate in glob("{0}/intermediates/{0}*".format(fam_name)):
        tar_file.add(intermediate)
        os.remove(intermediate)
    tar_file.close()
    # Move analysis to output directory
    os.makedirs(out, exist_ok=True)
    out_final = "{}/{}".format(out, fam_name)
    shutil.move(fam_name, out_final)


if __name__ == "__main__":
    # Define arguments
    parser = argparse.ArgumentParser("CTDG annotation pipeline")
    parser.add_argument("--name_family", "-n",
                        action='store',
                        help="Name of the gene family")

    parser.add_argument("--ref_seq", "-f",
                        action='store',
                        help="reference gene family sequence",
                        default=None)

    parser.add_argument("--ref_pfam", "-p",
                        action="store",
                        help="Reference Pfam Accession",
                        default=None)

    parser.add_argument("--hmmer_samples", "-S",
                        action="store",
                        type=int,
                        help="Number of samples to build empirical distribution of duplicates")
    parser.add_argument("--sp", "-s",
                        action="append",
                        default=[],
                        help="Species to run the analysis on (must be written in quotes)")
    parser.add_argument("--out_dir", "-o",
                        action="store",
                        default='CTDG_out',
                        help="Output directory")
    parser.add_argument("--cpu", "-c",
                        action="store",
                        type=int,
                        default=1,
                        help="CPU to be used")
    parser.add_argument("--db", "-d",
                        action="store",
                        help="Directory with gene annotation and HMMer database")

    parser.add_argument("--dir", "-D",
                        action="store",
                        default=None,
                        help="run analyses with all the sequences in a directory")
    # Check if hmmerscan is installed. Since this is not required for defining the analysis, it is executed before
    # the class definition
    hmmer_path = shutil.which('hmmscan')
    assert bool(hmmer_path), "HMMerScan is not installed or it is not in your PATH"

    # Parse arguments
    args = parser.parse_args()
    ctdg_config = CtdgConfig(out_dir=args.out_dir, db=args.db,
                             hmmer_samples=args.hmmer_samples, sp=args.sp)
    ctdg_config.check_arguments()
    ctdg_config.init_tables()
    # Initialize analysis
    analysis = CtdgRun(args=args, ctdg=ctdg_config)
    # Remove directory structure if partial results are found
    analysis.remove_if_unfinished()
    # Create directory structure
    analysis.create_folders()
    # Run HMMscan and Meanshift step
    hmmers = analysis.hmmscan()
    if hmmers:
        analysis.selected = parallel_meanshift(hmmers, args.cpu)
        # Organize summary table for each cluster candidate
        # pre_numbers = analysis.numbers.loc[:, ]
        # Run sampling
        numbers = sample_table(analysis)
        # Clean tables and save CSVs
        if numbers:
            clean_p95(numbers)
    else:
        print("No clusters were found")
    save_results(analysis)
    print("DONE")
