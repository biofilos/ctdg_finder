# include("ctdg_modules.jl")
# using CTDG
using DataFrames
using ZipFile
using ArgParse
using Base.Test
using JSON
using DataFrames
using Bio.Seq
using PyCall
using Formatting

@pyimport sklearn.cluster as cl

"""
Parse script options
"""
function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--name_family", "-n"
      help = "Name of gene family"
      required = false
    "--ref_seq", "-f"
      help = "Reference sequence (query)"
      required = false
      default = "None"
    "--hmm_ref", "-p"
      help="query HMM. ALL will run the analysis for all hmms in db"
      required= false
      arg_type = String
      default = "None"
    "--hmmer_samples", "-S"
      help = "Number of samples to build the empirical distribution"
      arg_type = Int
      default = 1000
    "--sp", "-s"
      help = "Restrict the analysis to a species set (optional)"
      nargs = '*'
      action = :append_arg
    "--out_dir", "-o"
      help = "Output directory"
      default = "CTDG_out"
    "--cpu", "-c"
      help = "Number of CPUs to use"
      arg_type = Int
      default = 1
    "--db", "-d"
      help = "Database to use"
      required = true
  end
  return parse_args(s)
end

function check_db(dir)
  int2str = x -> string(x)
  @assert ispath(dir) "Directory with CTDG database does not exist"
  db_files = readdir(dir)
  for i in ["pfam","chromosomes.csv","genes_parsed.csv", "hmmer"]
    @assert i in db_files "$i does not exist in $dir"
  end
  genomes = readtable(dir * "/chromosomes.csv")
  genes = readtable(dir * "/genes_parsed.csv")
  pfam_dict = JSON.parsefile(dir * "/hmmer/pfam.json")
  # Only use columns of interest
  genomes = genomes[:, [:species, :chromosome, :length, :bandwidth]]
  genes = genes[:, [:acc, :chromosome, :species, :symbol,
		    :start, :_end, :strand, :length]]
  genes[:, :strand] = map(int2str, genes[:, :strand])
  genes[:, :chromosome] = map(int2str, genes[:, :chromosome])
  return (genes, genomes, pfam_dict)
end

"""
Test if `HMMer` exists in PATH
"""
function test_hmmer()
  # return (./hmmsearch, true)
  hmmer_exist = false
  hmmer_path = "None"
  for path = split(ENV["PATH"],":")
    if ispath(path * "/hmmsearch") && ispath(path * "/hmmscan")
      hmmer_exist = true
      hmmer_path = path
      break
    end
  end
  return (hmmer_path, hmmer_exist)
end

"""
Remove unfinished analysis
"""
function remove_if_incomplete(name_family)
  if ispath(name_family)
    println("Partial results found for $name_family. Removing them")
    Base.rm(name_family, recursive=true)
  end
end

"""
Check that the selected species (if any) exist in the database
"""
function check_sp(args)
  if ! isempty(args["sp"])
    valid_species = unique(genomes[:species])
    for spp=args["sp"]
      @assert spp in valid_species "$spp is not a valid species"
    end
    println("Running the analysis with the species:")
    println(join(args["sp"], "\n"))
  end
end

"""
Create analysis directory structure
"""
function create_folders(name_family)
  for folder=["/intermediates", "/report"]
    mkpath(name_family * folder)
  end
end

"""
Saves the names of each pfam accession found in a hmmer search in a file
"""
function save_pfams(name_family, hmmer_out)
  println("processing Pfam output")
  pfams = []
  fileO = open(format("{1}/intermediates/{1}.pfams", name_family), "w")
  for line in eachline(open(hmmer_out))
    if line[1] != '#'
      data_line = [x for x in split(line, ' ') if length(x) > 0]
      pfam = data_line[2]
      push!(pfams, pfam)
      write(fileO, pfam*"\n")
    elseif contains(line, "[ok]")
      close(fileO)
      return pfams
  end
  end
end

"""
Execute HMMscan search
"""
function hmmscan(hmmer_path; cpus=1, ref="None",
                db="", evalue=1, pfam="None", genes=DataFrame(),
                pfam_dict=Dict(), name_family="")
  db_hmm = db * "/../pfam/pfam.hmm"
  all_seqs = db * "/seqs_hmmer.fa"
  if ref != "None"
    hmmer_scan = hmmer_path* "/hmmscan"
    out_file = format("{1}/intermediates/{1}.pre_hmmer", name_family)
    dom_out = out_file * "_dom"
    long_out = out_file * "_long"
    println("Running HMMscan with E-value = $evalue")
    run(`$hmmer_scan -o $long_out --tblout $out_file --domtblout $dom_out
    -E $evalue --cpu $cpus $db_hmm $ref`)
    # Save Pfam accessions on a file
    pfams = save_pfams(name_family, out_file)
  elseif pfam == "ALL":
    acc_file = db * "/../pfam/pfam.accs"
    pfams = [chomp(x) for x in readlines(acc_file)]
  else
    pfams = [pfam]
  end
  genes_dict = Dict()
  for pfam = pfams
    gene_list = pfam_dict[pfam]
    genes_slice = genes[findin(genes[:acc], gene_list), :]
    genes_dict[pfam] = genes_slice
  end
  return genes_dict
end

"""
Run the meanshift algorithm (prom Python scikit-learn) and annotate
relevant cluster candidates
"""
function MeanShift(hmmer_table, genomes, all_genes, pfam, name_family)
  hmmer_table[:, :family] = name_family
  hmmer_table[:, :pfam] = pfam
  hmmer_table[:, :cluster] = "0"
  hmmer_table[:, :order] = 0
  sort!(hmmer_table, cols=[:species, :chromosome, :start])
  for ms_sp_table=groupby(hmmer_table, [:species, :chromosome])
    sp_mean_shift = ms_sp_table[1, :species]
    chrom_mean_shift = ms_sp_table[1, :chromosome]

    if nrow(ms_sp_table) > 1
      bandwidth = genomes[(genomes[:species].==sp_mean_shift)&
                   (genomes[:chromosome].==chrom_mean_shift), :bandwidth][1]
      gene_starts = ms_sp_table[:, :start]
      gene_coords = [(x, y) for (x, y) in zip(gene_starts,
                                              zeros(Int, length(gene_starts)))]

      ms = cl.MeanShift(bandwidth=bandwidth)
      group_labels = ms[:fit_predict](gene_coords)

      ms_sp_table[:, :cluster] = map((x,y) -> "$(pfam)_" * string(x) * "_" * string(y),
      ms_sp_table[:chromosome], group_labels)
      for group in groupby(ms_sp_table, :cluster)
        if nrow(group) < 2
          group[:, :cluster] = "na_ms"
        else
          group[:, :order] = range(1, nrow(group))
        end
      end
    end
  end
  return hmmer_table
end


"""
Summarize information for each cluster
"""
function cluster_numbers(hmmer_table, all_genes)
  groups = groupby(hmmer_table, [:species, :chromosome, :cluster])
  n_groups = length(groups)
  # Initialize arrays
  sp_arr = Array(String, n_groups)
  chrom_arr = Array(String, n_groups)
  clust_arr = Array(String, n_groups)
  start_arr = Array(Int, n_groups)
  end_arr = Array(Int, n_groups)
  dups_arr = Array(Int,n_groups)
  total_arr = Array(Int, n_groups)
  prop_arr = Array(Float64, n_groups)
  # annotate each
  for (ix,group)=enumerate(groups)
    n_clustered = nrow(group)
    sp = group[1, :species]
    chrom = group[1, :chromosome]
    cluster = group[1, :cluster]
    start = minimum(group[:start])
    finish = maximum(group[:_end])
    all_cluster_tab = all_genes[(all_genes[:species] .== sp)&
    (all_genes[:chromosome] .== chrom)&
    (all_genes[:start] .>= start)&
    (all_genes[:_end] .<= finish), :]
    # not_clustered = length(setdiff(Set(all_cluster_tab[:acc]),
    #                                Set(group[:prot_acc])))
    total_genes = nrow(all_cluster_tab)
    prop_clustered = n_clustered / total_genes

    # Assign to the corresponding arrays
    sp_arr[ix] = sp
    chrom_arr[ix] = chrom
    clust_arr[ix] = cluster
    start_arr[ix] = start
    end_arr[ix] = finish
    dups_arr[ix] = n_clustered
    total_arr[ix] = total_genes
    prop_arr[ix] = prop_clustered
  end
  numbers = DataFrame(species=sp_arr,
  chromosome=chrom_arr,
  cluster=clust_arr,
  start=start_arr,
  _end=end_arr,
  duplicates=dups_arr,
  total_genes=total_arr,
  proportion_duplicates=prop_arr)
  numbers[:, :perc95_gw] = 0.0
  return numbers[numbers[:duplicates] .> 1, :]
end


"""
Runs Meanshift on each pfam table, and return a merged annotated table
with cluster candidates
"""
function run_ms(hmm_out, genomes, genes, name_family)
  ms_list, numbers_big = [], []
  for (pf_id, hmm_df) = hmm_out
    ms_table = MeanShift(hmm_out[pf_id], genomes, genes, pf_id, name_family)
    ms_numbers = cluster_numbers(ms_table, genes)
    push!(numbers_big, ms_numbers)
    push!(ms_list, ms_table)
  end
  return (ms_list, numbers_big)
end

"""
Given a list of accessions, return the maximum number of duplicates with a
specified E-value that are in that list (which correspond to all the genes
in a region
`acc_list`: all the genes in a region
`duplicates`: Dictionary of parsed blast in the chromosome
`evalue`: Target E-value
"""
f(x) = 1
@everywhere function grab_duplicates(acc_list, duplicates)
  duplicates_list = []
  duplicates_eval = Dict()
  # Use the accessions that have blast results
  acc_common = intersect(Set(acc_list),Set(keys(duplicates)))
  hit_counts = Dict()
  # Go through the genes in the region
  for acc in acc_common
    # Start hits counter for each accession
    duplicates_eval = 0
    # Each gene can have blast hits on several other
    # genes, go through them
    for hit in duplicates[acc]
      # If the gene has a blast hit to another gene in the same
      # region, count it.
      if hit  in acc_list
	       duplicates_eval += 1
      end
    end
    # Include the gene only if t
    if duplicates_eval > 1
      hit_counts[acc] = duplicates_eval
    end
  end
  if length(hit_counts) > 0
    return maximum(values(hit_counts))
  else
    return 0
  end
end


"""
Calculate the max-duplicates from one random sample, given a species table
`sp`: species name
`all_genes`: gene annotation
`chrom_table`: table containing chromosome lengths
`region_len`: length of the region being sampled
`evalue`: evalue threshold
"""
function sample_sp(numbers_row, all_genes, chrom_table, args)
  sp = numbers_row[:species]
  region_len = numbers_row[:_end] - numbers_row[:start]
  sp_table = all_genes[all_genes[:species] .== sp, :]
  chroms = chrom_table[(chrom_table[:species] .== sp)&
		       # Only include chromosomes that are
		       # longer to the region being sampled
		       (chrom_table[:length] .> region_len),
		       [:chromosome, :length]]
  cluster = numbers_row[:cluster]
  max_values = @sync @parallel vcat for i=1:args["hmmer_samples"]

    # Get random chromosome
    sample_up = 0
    if length(chroms[:chromosome]) == 0
      println(format("sp table: {}|region length: {}|species: {}", nrow(sp_table), region_len, sp))
    end
    sample_chrom = sample(chroms[:chromosome])

    # Calculate the up limit for sampling chromosome regions
    sample_space = chroms[chroms[:chromosome] .== sample_chrom, :length] - region_len
    # Calculate sample
    sample_up = sample(range(1, sample_space[1]))
    sample_start, sample_end = (sample_up, sample_up + region_len)
    # Get accessions for the genes in the sample region
    sample_accs = sp_table[(sp_table[:chromosome] .== sample_chrom)&
			   (sp_table[:start] .< sample_end)&
			   (sp_table[:_end] .> sample_start),:acc]
    json_file = readstring(open("$(args["db"])/hmmers/$(sp)_$(sample_chrom).json"))
    blast_dic = JSON.parse(json_file)
    max_dups = grab_duplicates(sample_accs, blast_dic)
    outfile = "$(args["name_family"])/intermediates/samples.coords"
    open(outfile, "a") do fileO
      log = "$(sp),$(cluster),$(sample_chrom),$(sample_start),$(sample_end),$(max_dups)\n"
      write(fileO, log)
    end
    #println(max_dups)
    max_dups
  end
  return percentile(max_values, 95)
end

"""
Calculate the maximum number of duplicates in the samples
specified by `blast_samples` and return the table with the
a column `perc95_gw`
"""
function sample_table(numbers, all_genes, genomes, name_family,  args)
  col_lengths = ""
  for (ix,col) = enumerate([:species, :cluster])
    spaces = [length(x) for x in numbers[:,col]]
    if length(spaces) > 0
      max_space = maximum(spaces)
    else
      max_space = 0
    end
    col_lengths *= "{:$(max_space + 2)s} "
  end
  #col_lengths *= "{:>12s} {:15s}"
  str_fmt = FormatExpr(col_lengths * "{:>12s}  {:15s}")
  numbers = numbers[numbers[:cluster] .!= "na_ms", :]
  if nrow(numbers) >0
    println("$(nrow(numbers)) cluster candidates")
    printfmtln(str_fmt, "Species", "Cluster", "Duplicates", "Percentile 95")
  end
  str_fmt = FormatExpr(col_lengths * "{:>12d}  {:>13.3f}")
  for row = eachrow(numbers)
    if row[:cluster] == "na_ms"
      p95 = 0
      color = :gray
    else
      p95 = sample_sp(row, all_genes, genomes, args)
      if row[:duplicates] >= p95
	       color = :green
      else
	       color = :red
      end
    end
    row[:perc95_gw] = p95
    msg = format(str_fmt, row[:species], row[:cluster],
		 row[:duplicates], row[:perc95_gw])
    print_with_color(color, msg * "\n")
  end
  clean_number = numbers[numbers[:duplicates] .>= numbers[:perc95_gw], :]
  clusters = unique(clean_number[:cluster])
  deleteat!(clusters, findin(clusters ,["na_ms","0","na_95"]))
  clean_numbers = clean_number[findin(clean_number[:cluster], clusters), :]
  for_removal = numbers[numbers[:duplicates] .< numbers[:perc95_gw], :]
  return numbers, for_removal, clean_numbers
end

"""
Clean intermediates
"""
function clean_intermediates(args)
  intermediates = readdir("$(args["name_family"])/intermediates/")
  w = ZipFile.Writer("$(args["name_family"])/intermediates/intermediates.zip")
  for file=intermediates
    f = ZipFile.addfile(w, file)
    file = "$(args["name_family"])/intermediates/"* file
    write(f, read(open(file)))
    Base.rm(file)
  end
  close(w)
end

"""
Move finished analysis to output directory
"""
function move_finished(args)
  if ! ispath(args["out_dir"])
    mkdir(args["out_dir"])
  end
  mv(args["name_family"], "$(args["out_dir"])/$(args["name_family"])")
end




hmm_path = "/usr/local/bin"

args = parse_commandline()
remove_if_incomplete(args["name_family"])
create_folders(args["name_family"])
genes, genomes, pfam_d = check_db(args["db"])

hmm_out = hmmscan(hmm_path,
cpus=args["cpu"], ref=args["ref_seq"], pfam=args["hmm_ref"],
db=args["db"], evalue=0.0001, genes=genes,
pfam_dict=pfam_d, name_family=args["name_family"])

ms_tables, numbers_95 = run_ms(hmm_out, genomes, genes, args["name_family"])

hmmer_final = DataFrame()
number_final = DataFrame()
clean_numbers_final = DataFrame()
for ix = 1:length(ms_tables)
  number = numbers_95[ix]
  ms_table = ms_tables[ix]
  dirty_number, for_removal, clean_numbers = sample_table(number, genes, genomes,
                                                 args["name_family"], args)
  ms_table[:order] = [string(x) for x in ms_table[:order]]
  for (sp,clu) = zip(for_removal[:species],for_removal[:cluster])
  ms_table[(ms_table[:species] .== sp)&
           (ms_table[:cluster] .== clu), :cluster] = "na_95"
  end
  # Reannotate the order for excluded clusters
  # This was done here because I did'nt want to change the type of
  # the column but at the very end
  if nrow(hmmer_final) == 0
    hmmer_final = ms_table
  else
    append!(hmmer_final, ms_table)
  end
  if nrow(number_final) == 0
    number_final = dirty_number
  else
    append!(number_final, dirty_number)
  end
  if nrow(clean_numbers_final) == 0
    clean_numbers_final = clean_numbers
  else
    append!(clean_numbers_final, clean_numbers)
  end
end
# Print numbers files
tab_file_name = format("{1}/report/{1}_numbers.csv", args["name_family"])
writetable(tab_file_name, number_final)
clean_num_file = replace(tab_file_name, "numbers", "numbers_clean")
writetable(clean_num_file, clean_numbers_final)

 hmmer_final[hmmer_final[:cluster] .== "na_95", :order] = "na_95"
 hmmer_final[hmmer_final[:cluster] .== "na_ms", :order] = "na_ms"
 gene_file = format("{1}/report/{1}_genes.csv",args["name_family"])
 writetable(gene_file, hmmer_final)
 # Filter out non-clustered genes
 gene_clean_file = replace(gene_file, "genes","genes_clean")
 # Remove non-clustered order values
 orders = unique(hmmer_final[:cluster])
 deleteat!(orders, findin(orders ,["na_ms","0","na_95"]))
 # Filter out non-clustered genes
 genes_clean = hmmer_final[findin(hmmer_final[:cluster], orders), :]
 writetable(gene_clean_file, genes_clean)
 # Compress and save intermediate files
 clean_intermediates(args)
 # Move analysis files to output directory
 move_finished(args)
