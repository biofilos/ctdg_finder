using DataFrames
using Formatting
using FastaIO

s_genes = readtable(ARGS[1])
seq_db = Dict(split(name,'|')[2] => seq for (name, seq) in FastaReader(ARGS[2]))


function get_seq(record)
  seq = seq_db[record[:acc]]
  name = format(">{}|{}|{}|{}|{}|{}\n{}\n", record[:species],
                                      record[:chromosome],
                                      record[:acc],
                                      record[:symbol],
                                      record[:cluster],
                                      record[:order],
                                      seq)
  return (name, record[:acc])
end
dones = []
open(ARGS[3], "w") do fasta_out
  for record = eachrow(s_genes)
    fasta, name = get_seq(record)

    if ! (name in dones)
      write(fasta_out, fasta)
      push!(dones, name)
    end
  end
end
