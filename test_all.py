import unittest
from ctdg_finder import (load_config, get_genes,
                         group_chrom, get_chrom_homologs,
                         get_chrom_info)
from HTSeq import GenomicFeature
from io import StringIO

class TestConfig(unittest.TestCase):
    config = load_config("./config.json")
    def test_load_config(self):
        self.assertIsInstance(self.config, dict)

class TestGenes(unittest.TestCase):
    def run_gene(self, gff, feature):
        return [x for x in get_genes(gff, feature)]
    
    base_gene = ["chr19","Annotation","gene",
                 "1", "100", ".","+",".",
                 "ID=gene;families=fam1,fam2"
                ]
    bad_gene = base_gene.copy()
    bad_gene[2] = "bad_feature"
    good_str = "\t".join(base_gene)
    bad_str = "\t".join(bad_gene)
    # Include three records, two well built, and one
    # to be removed
    good_good_bad = [good_str, bad_str, good_str]
    g_g_b_str = "\n".join(good_good_bad)
    good_gene = StringIO(good_str)
    bad_gene = StringIO(bad_str)
    g_g_b_gene = StringIO(g_g_b_str)
    # Include three genes, two with the correct feature
    # but one of them without families annotation
    g_g_bad_ann_b = [good_str,
                     good_str.replace("fam1,fam2","None"),
                     bad_str]
    g_g_bad_ann_b_str = "\n".join(g_g_bad_ann_b)
    g_g_bad_ann_b_gene = StringIO(g_g_bad_ann_b_str)

    
    feature = "gene"
    def test_gets_good_gene(self):
        
        gene = self.run_gene(self.good_gene, self.feature)[0]
        self.assertIsInstance(gene, GenomicFeature)

    def test_bad_feature(self):
        gene = self.run_gene(self.bad_gene, self.feature)
        n_genes = len(gene)
        self.assertEqual(n_genes, 0)
        
    def test_filter_bad_out(self):
        genes = self.run_gene(self.g_g_b_gene, self.feature)
        n_genes = len(genes)
        self.assertEqual(n_genes, 2)
    
    def test_bad_family(self):
        genes = self.run_gene(self.g_g_bad_ann_b_gene, self.feature)
        n_genes = len(genes)
        self.assertEqual(n_genes, 1)
        
    def test_no_fam(self):
        # Include three genes, one of them with no families
        # in its features
        g_g_nofam_bad = [
        self.good_str,
        self.good_str.replace(";families=fam1,fam2", ""),
        self.bad_str
        ]
        g_g_nofam_bad_str = "\n".join(g_g_nofam_bad)
        g_g_nofam_bad_gene = StringIO(g_g_nofam_bad_str)
        genes_no_fam = self.run_gene(g_g_nofam_bad_gene, self.feature)
        n_genes_no_fam = len(genes_no_fam)
        expected_genes = 1
        self.assertEqual(n_genes_no_fam, expected_genes)
    
    def test_bad_feature(self):
        g_g_nofam_bad = [
            self.good_str,
            self.good_str.replace(";families=fam1,fam2", ""),
            self.bad_str
            ]
        g_g_nofam_bad_str = "\n".join(g_g_nofam_bad)
        g_g_nofam_bad_gene = StringIO(g_g_nofam_bad_str)
        """
        Test if changing the feature works.
        """
        genes = self.run_gene(g_g_nofam_bad_gene, "bad_feature")
        n_genes = len(genes)
        self.assertEqual(n_genes, 1)
        
        
class TestGroupChroms(unittest.TestCase):
    feature = "gene"

    def test_right_numbers(self):
        base_gene1 = ["chr19","Annotation","gene",
                "1", "100", ".","+",".",
                "ID=gene;families=fam1,fam2"
            ]
        base_gene2 = ["chr1","Annotation","gene",
                    "1", "100", ".","+",".",
                    "ID=gene;families=fam1,fam2"
                    ]
        gene1_str = "\t".join(base_gene1)
        gene2_str = "\t".join(base_gene2)
        three_good = [gene1_str, gene2_str, gene1_str]
        three_good_str = "\n".join(three_good)
        three_good_genes = StringIO(three_good_str)
        gene_dict = group_chrom(three_good_genes, self.feature)
        n_chroms = sorted(list(gene_dict.keys()))
        expected_chroms = sorted(["chr1", "chr19"])
        self.assertListEqual(n_chroms, expected_chroms)
    
    def test_right_genes_per_chrom(self):
        base_gene1 = ["chr19","Annotation","gene",
                "1", "100", ".","+",".",
                "ID=gene;families=fam1,fam2"
            ]
        base_gene2 = ["chr1","Annotation","gene",
                    "1", "100", ".","+",".",
                    "ID=gene;families=fam1,fam2"
                    ]
        gene1_str = "\t".join(base_gene1)
        gene2_str = "\t".join(base_gene2)
        three_good = [gene1_str, gene2_str, gene1_str]
        three_good_str = "\n".join(three_good)
        three_good_genes = StringIO(three_good_str)
        gene_dict = group_chrom(three_good_genes, self.feature)
        chrom_n = {x: len(y) for x, y in gene_dict.items()}
        expected_dict = dict(chr19=2,chr1=1)
        self.assertDictEqual(chrom_n, expected_dict)
        

class TestChromHomologs(unittest.TestCase):
    def test_good_groups(self):
        feature = "gene"
        base_gene1 = ["chr19","Annotation","gene",
                    "1", "100", ".","+",".",
                    "ID=gene;families=fam1,fam2"
                    ]
        gene1_str = "\t".join(base_gene1)
        gene2_str = gene1_str.replace("fam1,fam2", "fam1")
        gene3_str = gene1_str.replace("fam1,fam2", "fam3")
        gene4_str = gene1_str.replace("chr19","chr20")


        genes_file = StringIO("\n".join([gene1_str, gene2_str,
        gene3_str, gene4_str]))
        chrom_homologs = get_chrom_homologs(genes_file, feature)
        chr19_fam1_has_two = len(chrom_homologs["chr19"]["fam1"]) == 2
        chr19_fam2_has_one = len(chrom_homologs["chr19"]["fam2"]) == 1
        chr19_fam3_has_one = len(chrom_homologs["chr19"]["fam3"]) == 1
        good_fams_19 = chr19_fam1_has_two and chr19_fam2_has_one and chr19_fam3_has_one
        # There is only one gene in chr20, which belongs to two families
        chr20_fam1_has_one = len(chrom_homologs["chr20"]["fam1"]) == 1
        chr20_fam2_has_one = len(chrom_homologs["chr20"]["fam2"]) == 1
        good_fams_20 = chr20_fam2_has_one and chr20_fam1_has_one
        self.assertTrue(good_fams_19 and good_fams_20)



if __name__ == "__main__":
    unittest.main()
    