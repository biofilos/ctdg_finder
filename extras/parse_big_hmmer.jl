import JSON

pfam_in = ARGS[1]
dict_out = ARGS[2]


pfam_db = Dict()
for line = eachline(open(pfam_in))
  if line[1] != '#'
    data_line = [x for x in split(line, " ") if x != ""]
    seq  = split(data_line[1],'|')[2]
    pfam = data_line[4]
    if ! haskey(pfam_db, pfam)
      pfam_db[pfam] = [seq]
    else
      push!(pfam_db[pfam], seq)
    end
  end
end

stringdata = JSON.json(pfam_db)

# write the file with the stringdata variable information
open(dict_out, "w") do f
        write(f, stringdata)
     end
