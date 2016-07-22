from BCBio import GFF
# import pandas as pd
# from glob import glob
# from pprint import pprint
in_file = 'ref_GRCh38.p7_top_level.gff3'

limit_info = dict(gff_type=[('CDS',)])

# features_list = []

# for in_file in glob("*.gff3"):
in_handle = open(in_file)
# examiner = GFF.GFFExaminer()
# pprint(examiner.available_limits(in_handle))
# in_handle.close()
cds_file = open("genes_cds.csv", "w")
cds_file.write(','.join(['start', 'end', 'length', 'strand', 'acc', 'chromosome', 'species', 'symbol', 'name'])+'\n')
for rec in GFF.parse(in_handle, limit_info=limit_info, target_lines=1000):
    sp = in_file.replace("_scaffolds.gff3", '')
    chromosome_feat = rec.id
    for feat in rec.features:
        if feat.type == 'CDS':
            if 'gene' in feat.qualifiers:
                symbol_feat = feat.qualifiers['gene'][0]
            else:
                symbol_feat = feat.qualifiers['protein_id'][0]
            start_feat = feat.location.start.position
            end_feat = feat.location.end.position
            length_feat = abs(end_feat - start_feat)
            strand_feat = feat.location.strand
            acc_feat = feat.qualifiers['protein_id'][0]
            name = feat.qualifiers['product'][0]
            cds_file.write(",".join([start_feat, end_feat, length_feat, strand_feat, acc_feat,
                                     chromosome_feat, sp, symbol_feat, name]) + '\n')
            # features_list.append([start_feat, end_feat, strand_feat, acc_feat,
            #                       chromosome_feat, sp, symbol_feat, name])
in_handle.close()
# feat_table = pd.DataFrame(features_list, columns=['start', 'end', 'strand', 'acc', 'chromosome',
#                                                   'species', 'symbol', 'name']).set_index('acc')
# feat_table['length'] = abs(feat_table['start'] - feat_table['end'])
