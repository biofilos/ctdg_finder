import argparse
import json
import os
import shutil
import tarfile
from concurrent import futures
from functools import partial
from glob import glob
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.cluster import MeanShift
pd.options.mode.chained_assignment = None

# Define class 

def set_db(args):
    arg_dic = dict(sp=args.sp,
                   samples=args.hmmer_samples,
                   db=args.db,
                   out=args.out_dir,
                   genes_f=args.db + "/genes_parsed.csv",
                   chrom_d_f=args.db + "/chromosomes.json",
                   name=args.name_family,
                   ref_seq=args.ref_seq,
                   ref_pfam=args.ref_pfam,
                   cpu=args.cpu,
                   chrom_only=args.chrom_only)
    return arg_dic


def check_arguments(args):
    """
    Perform assertion statements to check if all the files needed
    for the analysis are present
    :return: None
    """
    # Check that the directory with the database exists
    assert os.path.exists(args["db"]), "directory {} with database does not exist".format(args["db"])
    #db_path = '{}/pfam/pfam.hmm'.format(self.db)
    db_path = '{}/panther/panther.hmm'.format(args["db"])
    # Check that the blast database files exists
    hmmer_files = sum([os.path.exists(db_path + x) for x in ['', '.h3f', '.h3i', '.h3m', '.h3p']])
    assert hmmer_files == 5, "Some of the HMMer database files with prefix pfam do not exist"
    # Check that the chromosomes information table exists
    assert os.path.exists(args["chrom_d_f"]), "Chromosome information file" + \
                                             "{}/chromosomes.json does not exist".format(args["db"])
    # Check that the annotation table exists
    assert os.path.exists(args["genes_f"]), "Annotation file" + \
                                             "{}/genes_parsed.csv does not exist".format(args["db"])

    # Check that number of blast samples has been set
    assert args["samples"], "Number of Hmmer samples was not specified"
    # Make sure that the output directory is not the directory where the database is or the directory    # the CTDGFinder resides
    banned_dirs = [os.getcwd(), args["db"], '.']
    assert args["out"] not in banned_dirs, "output direcrory can't be the ctdgfinder directory" + \
                                                " or the directory where the ctdg database is stored"

def genes_chroms(genes_f, chrom_d_f):
    gene_tab = pd.read_csv(genes_f)
    gene_tab.set_index("acc", inplace=True)
    gene_tab.loc[:, "chromosome"] = gene_tab["chromosome"].astype(str)

    chrom_d = json.load(open(chrom_d_f))
    return (gene_tab, chrom_d)


def create_folders(name):
    """
    Create the directories where all the
    results will be saved
    """
    if not os.path.exists(name):
        os.makedirs(name)
        folders = ['intermediates', 'report']
        for sub_dir in folders:
            target = name + '/' + sub_dir
            if not os.path.exists(target):
                os.makedirs(target)
    print("Directory tree for {} was created.".format(name))

def remove_if_unfinished(out, name):
    """
    Remove unfinished analysis
    :return: None
    """
    # Do not run analysis if final output is present
    assert not os.path.exists("{}/{}".format(out, name)), "Results for {} exist.".format(name)

    if os.path.exists(name):
        print("Partial results for {} were found, removing them".format(name))
        shutil.rmtree(name)

def strict_hmmer(hmmer_out, threshold):
    """
    filter out hmmer with less mutual coverare specified by the threshold argument
    :param threshold: proportion of sequence length (in query and subject) to be used as filter
    :param hmmer_out: output (.dom) of hmmerscan
    :return: hmmer accessions that passed the threshold
    """
    pfam_filtered_fn = hmmer_out.replace("_dom", "_filtered_dom")
    pfam_filtered = open(pfam_filtered_fn, "w")
    pfams = []
    for line in open(hmmer_out):
        # Skip comment lines
        if not line.startswith("#"):
            # Process each line
            data_line = [x for x in line.split(' ') if x != ""]
            acc = data_line[0].split("|")[0]
            pfam = data_line[3]
            # Include entries with a minimum coverage
            # in both the hmm and the query of 30%
            hmm_start = int(data_line[15])
            hmm_end = int(data_line[16])
            hmm_aln_len = hmm_end - hmm_start
            hmm_len = int(data_line[5])
            hmm_coverage = hmm_aln_len / hmm_len

            q_len = int(data_line[2])
            q_end = int(data_line[20])
            q_start = int(data_line[19])
            q_aln_len = q_end - q_start
            query_coverage = q_aln_len / q_len
            if query_coverage >= threshold and hmm_coverage >= threshold:
                pfam_filtered.write(",".join(data_line))
                pfams.append(pfam)
    pfam_filtered.close()
    return pfams

def hmmscan(args, evalue=1e-3):
    name = args["name"]
    ref_pfam = args["ref_pfam"]
    ref_seq = args["ref_seq"]
    db = args["db"]
    cpu = args["cpu"]
    out_file = "{0}/intermediates/{0}.pre_hmmer".format(name)
    dom_out = out_file + "_dom"
    long_out = out_file + "_long"
    # pfam_hmm = self.db + "/pfam/pfam.hmm"
    pfam_hmm = db + "/panther/panther.hmm"
    if ref_seq and not ref_pfam:
        print("Running HMMscan")
        cmd = "hmmsearch -o {} --tblout {} --domtblout {} -E {} --cpu {} {} {}".format(long_out, out_file, dom_out, evalue, cpu, pfam_hmm, ref_seq)
        # Run HMMerScan
        os.system(cmd)

        # Parse the Pfam accessions found
        pfams = strict_hmmer(dom_out, 0.3)
    else:
        pfams = [ref_pfam]
    # Load Pfam Dictionary
    pfam_dict = json.load(open("{}/hmmer/pfam.json".format(db)))
    # Filter gene table with the genes that contain the Pfam domains found by HMMscan
    pfam_genes = []
    for pfam in pfams:
        if pfam in pfam_dict:
            pfam_genes += pfam_dict[pfam]
    # Remove duplicates
    pfam_selected = set(pfam_genes)
    num_candidates = len(pfam_selected)
    genes_f = open(args["genes_f"])
    columns = next(genes_f).strip().split(",")
    # Check that columns are in the right order
    msg = """
    Check that columns are in this order:
        acc, species, chromosome, start, end, strand
    """
    assert columns == ["acc", "species", "chromosome", "start", "end",
                       "strand"], msg
    genes_dict = dict()
    candidates = 0
    for gene in genes_f:
        acc, species, chromosome, start, end, strand = gene.split(",")
        bandwidth = chrom_d[species][chromosome][1]
        if acc in pfam_selected:
            # fields:
            # acc, start, end, cluster, order
            gene_coords = [acc, int(start), int(end), 0, 0]
            if species not in genes_dict:
                genes_dict[species] = {chromosome: [bandwidth,
                                                    [gene_coords]]}
            elif chromosome not in genes_dict[species]:
                    genes_dict[species][chromosome] = [bandwidth,
                                                       [gene_coords]]
            else:
                genes_dict[species][chromosome][1].append(gene_coords)
            candidates += 1
            pfam_selected.remove(acc)
        if candidates == num_candidates:
            break
    return genes_dict


def setup_meanshift(sp_dictionary, name):
    """
    Parses dictionary of the form 
    {sp[chrom]:[bandwidth,[acc,start,end,cluster,order]]
    and run MeanShift on each one 
    """
    for sp in sp_dictionary:
        for chrom in sp_dictionary[sp]:
            bw, genes = sp_dictionary[sp][chrom]
            # Sort genes by their start coordinate
            genes = sorted(genes, key=lambda x: x[1])
            # Process chromosomes with more than 1
            # gene from the target gene family
            if len(genes) > 1:
                ms_annotated = meanshift(genes, bw)
                counts = Counter(ms_annotated[1])
                # Flag cluster candidates for removal if they
                # only have one gene
                ms_remove = [x for x in counts if counts[x] == 1]
                order = 1
                past_cl = ""
                n_genes = len(ms_annotated[0])
                annotated = [0] * n_genes
                count = 0
                for gene, label in zip(*ms_annotated):
                    cl_name = "{}_{}_{}".format(name, chrom,
                                                str(label + 1))
                    # If new cluster, reset coordinates
                    if cl_name != past_cl:
                        order = 1
                        past_cl = cl_name
                    # Flag genes outside cluster candidates
                    # at the meanshift step
                    if label in ms_remove:
                        cl_name = "na_ms"
                        order = 0
                    gene[3] = cl_name
                    gene[4] = order
                    order += 1
                    annotated[count] = gene
                    count += 1
            else:
                annotated = genes
            sp_dictionary[sp][chrom] = annotated
    return sp_dictionary

def meanshift(gene_coords, bandwidth):
    """
    Run meanshift for a list of gene coordinates (start, end)
    """
    starts = [x[1] for x in gene_coords]
    gene_starts = np.array([(x, y) for x, y in zip(starts,
                                                   np.zeros(len(starts)))])
    mean_shift = MeanShift(bandwidth=bandwidth)
    try:
        mean_shift.fit(gene_starts)
    except ValueError:
        print("BW: {}".format(bandwidth))
    labels = mean_shift.labels_
    return (gene_coords, labels)


def build_numbers(ms_annotation):
    """
    Compile data for stats sampling (numbers table)
    """
    # Initialize numbers dictionary
    #numbers_dict = {sp: {} for sp in ms_annotation.keys()}
    numbers_lst = []
    for sp in ms_annotation:
        for chrom in ms_annotation[sp]:
            c_data = ms_annotation[sp][chrom]
            #if len(chrom_data) > 1:
            # Initialize data structure
            clusters = {x[3]: [] for x in c_data}
            if len(clusters) > 1:
                #if chrom not in numbers_dict[sp]:
                #    numbers_dict[sp][chrom] = {}
                for cluster in clusters:
                    min_coord = min([x[1] for x in c_data if cluster in x[3]])
                    max_coord = max([x[2] for x in c_data if cluster in x[3]])
                    num_dupli = len([x for x in c_data if cluster in x[3]])
                    numbers_lst.append([sp, chrom, cluster, min_coord,
                                                        max_coord,
                                                        num_dupli,
                                                        # Reserve for P95
                                                        0])
    numbers_lst.sort(key=lambda x: x[0])
    return numbers_lst





def build_df(numbers, genes, name):
    n_df = pd.DataFrame(numbers)
    # Linearize genes
    linear_genes = []
    for sp, chrom_data in genes.items():
        for chrom_name, genes_data in chrom_data.items():
            for gene in genes_data:
                linear_genes.append([sp, chrom_name] + gene)
    g_df = pd.DataFrame(linear_genes)

    # Set column names
    n_df.columns = ("species", "chromosome", "cluster", "start",
                    "end", "duplicates", "P95")

    g_df.columns = ("species", "chromosome", "acc", "start",
                    "end", "cluster", "order")
    # Get accessions that did not pass the percentile threshold
    for_removal = n_df["duplicates"] < n_df["P95"]
    na_95_tag = "na_95"
    for sp, ch, cl, *rst in n_df.loc[for_removal].set_index("species").to_records():
        g_df.loc[(g_df["species"] == sp) &
                 (g_df["cluster"] == cl), "cluster"] = na_95_tag
    n_df.loc[for_removal, "cluster"] = na_95_tag
    # Save tables
    discarded = ["0", "na_95", "na_ms", 0]

    n_file = "{0}/report/{0}_numbers.csv".format(name)
    g_file = n_file.replace("numbers", "genes")
    n_clean_file = n_file.replace("numbers", "numbers_clean")
    g_clean_file = n_clean_file.replace("numbers", "genes")

    n_df.to_csv(n_file)
    n_df.loc[~(n_df["cluster"].isin(discarded))].to_csv(n_clean_file)

    g_df.to_csv(g_file)
    g_df.loc[~(g_df["cluster"].isin(discarded))].to_csv(g_clean_file)



def get_p95(numbers, args, chrom_d, genes):
    cpu = args["cpu"]
    longest_sp = max([len(x[0]) for x in numbers])
    longest_cl = max([len(x[2]) for x in numbers])
    len_data = sum([1 for x in numbers if x[2] != "na_ms"])
    longest_sp += 2
    longest_cl += 2
    print("Sampling {} cluster candidates".format(len_data))
    print_msg = "{:<{sp}} {:<{cluster}} {:<10} {}{}"
    print(print_msg.format("Species", "Cluster", "Duplicates", "P95", "",
                           sp=longest_sp, cluster=longest_cl))
    with futures.ProcessPoolExecutor(cpu) as pool:
        new_numbers = [None] * len(numbers)
        count = 0
        for chrom_pack in pool.map(part_sample_chrom, numbers):
            data = numbers[count]
            if not data[2] == "na_ms":
                p95 = np.percentile(chrom_pack, 95)
                data[-1] = p95.round(4)
                if data[-2] >= data[-1]:
                    msg = "*"
                else:
                    msg = ""
                print(print_msg.format(data[0], data[2], data[-2],
                                       data[-1], msg, sp=longest_sp,
                                       cluster=longest_cl))
            new_numbers[count] = data
            count += 1
    return new_numbers

def get_sample_genes(chrom, sp, length, genes_df, db):
    sample_len = chrom_d[sp][chrom][0]
    min_coord = 0
    max_coord = sample_len - length
    random_start = np.random.randint(min_coord, max_coord)
    random_end = random_start + length
    genes = set(genes_df.loc[(genes_df["species"] == sp) &
                             (genes_df["chromosome"] == chrom) &
                             (genes_df["start"] < random_end) &
                             (genes_df["end"] > random_start)].index.values)
    mat_file = "{}/matrices/{}_{}.csv".format(db, sp, chrom)
    #return (pd.read_csv(mat_file,index_col=0), genes)
    return (mat_file, genes)

def get_max_dups(mat_genes):
    mat_file, genes = mat_genes
    #genes = genes_region(genes_df, sp, chrom, random_start, random_end)
    if len(genes) <= 1:
        max_dups = 0
    else:
        chrom_df = pd.read_csv(mat_file, index_col=0)
        in_both = genes.intersection(set(chrom_df.index.values))
        max_dups = chrom_df.loc[in_both, in_both].sum().max()
        if max_dups > 0:
            max_dups += 1
        if max_dups is np.nan:
            max_dups = 0
    return max_dups

def sample_chromosomes(cluster_data, samples, only_chrom, chrom_d, genes_df,
                       db, cpu, percent=95):
    sp, chromosome, cl_name, cl_start, cl_end, dups, p = cluster_data
    if cl_name == "na_ms":
        return cluster_data
    length = cl_end - cl_start
    # Only include chromosomes at least as long as the length of the sample
    if only_chrom:
        chromosomes = [chromosome] * samples
    else:
        chromosomes = tuple([x for x in chrom_d[sp] if
                             chrom_d[sp][x][0] >= length])

    chroms_samples = np.random.choice(chromosomes, samples)
    chrom_pack = [None] * samples
    chrom_pack = [get_max_dups(get_sample_genes(x, sp, length, genes_df, db)) for x in
                  chroms_samples]
    return chrom_pack



def genes_region(genes_df, sp, chrom, start, end):
    genes = set(genes_df.loc[(genes_df["species"] == sp) &
                             (genes_df["chromosome"] == chrom) &
                             (genes_df["start"] < end) &
                             (genes_df["end"] > start)].index.values)
    return genes



def save_results(args):
    # Compress intermediates
    out = args["out"]
    fam_name = args["name"]
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


       # Load the genes that have belong to the identified gene families
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
parser.add_argument("--chrom_only", "-C",
                    action="store_true",
                    help="Use samples only from the chromosome where the \
                    cluster candidate is found")

# Check if hmmerscan is installed. Since this is not required for defining the analysis, it is executed before
# the class definition
hmmer_path = shutil.which('hmmscan')
assert bool(hmmer_path), "HMMerScan is not installed or it is not in your PATH"

# Parse arguments
args = parser.parse_args()
data = set_db(args)
check_arguments(data)
remove_if_unfinished(data["out"], data["name"])
# Set reference data
genes, chrom_d = genes_chroms(data["genes_f"], data["chrom_d_f"])
create_folders(data["name"])
# Identify gene families
hmmers = hmmscan(data)
# run MeanShift
ms_annotation = setup_meanshift(hmmers, data["name"])
# Summarize cluster candidates
numbers = build_numbers(ms_annotation)
if numbers:
    part_sample_chrom = partial(sample_chromosomes, samples=data["samples"],
                                only_chrom=data["chrom_only"],
                                genes_df=genes,
                                db=data["db"],
                                cpu=data["cpu"],
                                chrom_d=chrom_d)
    new_numbers = get_p95(numbers, data, chrom_d, genes)
    build_df(new_numbers, ms_annotation, data["name"])
save_results(data)
print("DONE")
