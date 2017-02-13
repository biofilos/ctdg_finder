using DataFrames


"""
Given a database directory, calculate the bandwidth parameter for
MeanShift and include it in the chromosomes.csv file
"""
db_dir = ARGS[1]

function intergenic(chrom_table)
  # If there is only one gene in the chromosome, the "intergenic distance"
  # will be the end coordinate of the gene
  if nrow(chrom_table) == 1
    bandwidth = chrom_table[:_end][1]
    # If there are two genes, ID will be the end coordinate of the last
  elseif nrow(chrom_table) == 2
    bandwidth = abs(chrom_table[2, :start] - chrom_table[1, :_end])
  else
    pre_distances = chrom_table[2:end, :start] - chrom_table[1:end-1, :_end]
    gene_distances = pre_distances[pre_distances .> 0]
    if length(gene_distances) == 0
      bandwidth = sort(chrom_table[:_end])[end]
    elseif length(gene_distances) == 1
      bandwidth = gene_distances[1]
    else
      bandwidth = mean(gene_distances) + std(gene_distances)
    end
  end
  return bandwidth
end

genes = readtable(db_dir * "/genes_parsed.csv")
chromosomes = readtable(db_dir * "/chromosomes.csv")
sort!(genes,cols=[:species, :chromosome, :start])


intergenics = by(genes,[:species, :chromosome]) do df
  DataFrame(bandwidth=intergenic(df))
end
joined = join(chromosomes, intergenics, on=[:species,
                                            :chromosome])[:,[:species,
                                                             :chromosome,
                                                             :length,
                                                             :bandwidth]]
writetable(db_dir * "/chromosomes.csv", joined)
