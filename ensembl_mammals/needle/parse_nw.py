import sys

input_file = sys.argv[1]

def clean_line(line):
    """
    Takes a line from needleall (from EMBOSS) that contains a similarity or identity score,
    extracts the relevant number and convert it to a decimal fraction
    :param line: line containing identity or similarity
    :return: decimal fraction of identity or similarity
    """
    cleaned = line.split("(")[-1].strip(' ').replace("%)", '')
    clean_num = float(cleaned) / 100
    return clean_num

parsed_dict = {}

for line in open(input_file).readlines():
    if line.startswith("#"):
        if "# 1:" in line:
            query = line.split("|")[2]
        if "# 2:" in line:
            subject = line.split("|")[2]
        if "# Identity:" in line:
            identity = clean_line(line)
        if "# Similarity:" in line:
            similarity = clean_line(line)

        data_list = [subject, identity, similarity]
        if "Gaps" in line:
            if query not in parsed_dict:
                parsed_dict[query] = data_list
            else:
                parsed_dict[query].append(data_list)