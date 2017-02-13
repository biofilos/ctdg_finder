using DataFrames
using Formatting
using FastaIO

s_genes = readtable(ARGS[1])
seq_db = Dict(split(name,'|')[3] => seq for (name, seq) in FastaReader(ARGS[2]))
fasta_out =


function get_seq(record)
  seq = seq_db[record[:prot_acc]]
  name = format(">{}|{}|{}|{}|{}|{}\n{}\n", record[:species],
                                      record[:chromosome],
                                      record[:prot_acc],
                                      record[:symbol],
                                      record[:cluster],
                                      record[:order],
                                      seq)
  return name
end

open(ARGS[3], "w") do fasta_out
  for record = eachrow(s_genes)
    write(fasta_out, get_seq(record))
  end
end
