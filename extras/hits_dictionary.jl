using JSON
using DataFrames
using Formatting

db = ARGS[1]
genes_file = db * "/genes_parsed.csv"

pfams = JSON.parsefile(format("{}/hmmer/pfam.json", db))

gene_pfam = DataFrame(pfam=String[], acc=String[])

for (pfam, genes) = pfams
  for gene = genes
    push!(gene_pfam, [pfam, gene])
  end
end

genes_df = readtable(genes_file)

genes_annotated = join(genes_df, gene_pfam, on=:acc, kind=:left)

for chrom_table = groupby(genes_annotated, [:species, :chromosome])
  species = chrom_table[:species][1]
  chrom = chrom_table[:chromosome][1]
  chrom_file = format("{}/hmmers/{}_{}.json", db, species, chrom)
  chrom_dict = Dict()
  for gene in unique(chrom_table[:acc])
    gene_df = chrom_table[chrom_table[:acc].==gene, :]
    gene_df = gene_df[~(isna(gene_df[:pfam])),:]
    pfams_accs = chrom_table[(findin(chrom_table[:pfam],gene_df[:pfam])),:]
    chrom_dict[gene] = unique(pfams_accs[:acc])
  end
  dict_str = JSON.json(chrom_dict)
  fileO = open(chrom_file, "w")
  write(fileO, dict_str)
  close(fileO)
end
