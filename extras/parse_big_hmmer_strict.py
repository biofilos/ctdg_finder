import json
from sys import argv
"""
Parse table output from HMMer and save it as a JSON dictionary
"""


def main(pfam_in, dict_out_f, filtered_out, return_dict=False):
    # Initialize pfam dictionary
    pfam_db = {}
    dict_out = open(dict_out_f, "w")

    pfam_filtered = open(filtered_out, "w")
    for line in open(pfam_in):
        # Skip comment lines
        if not line.startswith("#"):
            # Process each line
            data_line = [x for x in line.split(' ') if x != ""]
            acc = data_line[0]
            seq = acc
            pfam = data_line[3]
            # Include entries with a minimum coverage
            # in both the hmm and the query of 30%
            threshold = 0.3
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
                if pfam not in pfam_db:
                    pfam_db[pfam] = [seq]
                else:
                    pfam_db[pfam].append(seq)

    # Save to file
    dict_out.write(json.dumps(pfam_db))
    dict_out.close()
    pfam_filtered.close()
    if return_dict:
        return pfam_db

if __name__ == "__main__":
    # Pfam tbl path
    pfam_in = argv[1]
    # output path
    dict_out = argv[2]
    filtered_out = argv[3]
    main(pfam_in, dict_out, filtered_out)

