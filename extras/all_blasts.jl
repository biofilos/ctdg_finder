using DataFrames
using Formatting
using Gadfly
"""
Parse blast output using bitscore
"""
function parse_bl_bs(out_file)
  # Load blast output (assumes query and subject are the same species)
  blast = readtable(out_file,
                    header=false,
                    separator='\t',
                    names=[:query, :subject, :qlen, :slen,
                           :qstart, :qend, :sstart, :send,
                           :length, :gaps, :gapopen, :evalue, :bitscore])
  # Extracts the bitscore of each sequence against itself
  reference = blast[blast[:query] .== blast[:subject], [:query, :bitscore]]
  bitscore_dict = Dict()
  for row in range(1,nrow(reference))
    data_row = reference[row, :]
    bitscore_dict[data_row[:query][1]] = data_row[:bitscore][1]
  end
    return (blast, bitscore_dict)
end
cd("/home/alekos/Documents/phd/vanderbilt/clusters/ctdg_db_process/blasts/")
blast, bitscore_dict = parse_bl_bs("chimp_chimp.blast")

blast_parsed = DataFrame(seq1=String[], seq2=String[], norm_bitscore=Float32[])
for seq1 in unique(blast[:, :query])
  # Get bitscore against itself for sequence 1
  seq1_bit = bitscore_dict[seq1]
  # Get the list of sequences that had a hit to seq1
  seq1_subjects = blast[blast[:query] .== seq1, :subject]
  for seq2 in unique(seq1_subjects)
    # Get bitscore against itself for sequence 2
    seq2_bit = bitscore_dict[seq2]
    # Get blist of sequences that had a hit with seq2
    seq2_subjects = blast[blast[:query] .== seq2, :subject]
    # Get bitscore of first sequence as query and seq2 as subject
    # and vice versa
    if seq1 != seq2 && seq1 in seq2_subjects
      hit12 = blast[(blast[:query].==seq1)&
                    (blast[:subject].==seq2),:bitscore][1]
      hit21 = blast[(blast[:query].==seq2)&
                    (blast[:subject].==seq1),:bitscore][1]

      normalized_bitscore = (hit12 + hit21) / (seq1_bit + seq2_bit)
      push!(blast_parsed, [seq1 seq2 normalized_bitscore])
      println(format("{} - {}: {}", seq1, seq2, normalized_bitscore))
    end
  end
end

blast_parsed = readtable("parsed_blast.csv")
# Include Evalue
rename!(blast_parsed, [:seq1,:seq2], [:query,:subject])
with_e = join(blast_parsed, blast,on=[:query,:subject],kind=:left)
with_e = with_e[:, [:query, :subject, :norm_bitscore, :evalue,:bitscore]]
with_e[:,:log_bs] = log2(with_e[:norm_bitscore])
