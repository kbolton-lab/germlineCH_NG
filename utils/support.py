# %%
import csv
import numpy as np
import os

# Amino Acid Code:
Amino_acid_abbreviations = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "*": "*",
}

Genetic_Code = {
    "TTT": "Phe",
    "TTC": "Phe",
    "TTA": "Leu",
    "TTG": "Leu",
    "TCT": "Ser",
    "TCC": "Ser",
    "TCA": "Ser",
    "TCG": "Ser",
    "TAT": "Tyr",
    "TAC": "Tyr",
    "TAA": "*",
    "TAG": "*",
    "TGT": "Cys",
    "TGC": "Cys",
    "TGA": "*",
    "TGG": "Trp",
    "CTT": "Leu",
    "CTC": "Leu",
    "CTA": "Leu",
    "CTG": "Leu",
    "CCT": "Pro",
    "CCC": "Pro",
    "CCA": "Pro",
    "CCG": "Pro",
    "CAT": "His",
    "CAC": "His",
    "CAA": "Gln",
    "CAG": "Gln",
    "CGT": "Arg",
    "CGC": "Arg",
    "CGA": "Arg",
    "CGG": "Arg",
    "ATT": "Ile",
    "ATC": "Ile",
    "ATA": "Ile",
    "ATG": "Met",
    "ACT": "Thr",
    "ACC": "Thr",
    "ACA": "Thr",
    "ACG": "Thr",
    "AAT": "Asn",
    "AAC": "Asn",
    "AAA": "Lys",
    "AAG": "Lys",
    "AGT": "Ser",
    "AGC": "Ser",
    "AGA": "Arg",
    "AGG": "Arg",
    "GTT": "Val",
    "GTC": "Val",
    "GTA": "Val",
    "GTG": "Val",
    "GCT": "Ala",
    "GCC": "Ala",
    "GCA": "Ala",
    "GCG": "Ala",
    "GAT": "Asp",
    "GAC": "Asp",
    "GAA": "Glu",
    "GAG": "Glu",
    "GGT": "Gly",
    "GGC": "Gly",
    "GGA": "Gly",
    "GGG": "Gly",
}

base_alternative = {"A": ["T", "C", "G"], "T": ["A", "C", "G"], "G": ["A", "C", "T"], "C": ["A", "G", "T"]}

Possible_AA_Changes = {}
Possible_Base_Changes = {}
Possible_Codon_Changes = {}
for k, v in Genetic_Code.items():
    codon = k
    base_1 = str(codon[0])
    base_2 = str(codon[1])
    base_3 = str(codon[2])
    new_codons = []
    amino_acid_changes = []
    base_changes = []
    if base_1 in base_alternative.keys():
        alternate_base = base_alternative[base_1]
        for i in alternate_base:
            base_changes.append(i)
            new_codon = i + base_2 + base_3
            new_codons.append(new_codon)
    if base_2 in base_alternative.keys():
        alternate_base = base_alternative[base_2]
        for i in alternate_base:
            base_changes.append(i)
            new_codon = base_1 + i + base_3
            new_codons.append(new_codon)
    if base_3 in base_alternative.keys():
        alternate_base = base_alternative[base_3]
        for i in alternate_base:
            base_changes.append(i)
            new_codon = base_1 + base_2 + i
            new_codons.append(new_codon)
    Possible_Codon_Changes[k] = (
        new_codons[0],
        new_codons[1],
        new_codons[2],
        new_codons[3],
        new_codons[4],
        new_codons[5],
        new_codons[6],
        new_codons[7],
        new_codons[8],
    )  # i.e. {'TTT': ['ATT', 'CTT', 'GTT'...etc]}
    Possible_Base_Changes[k] = (
        base_changes[0],
        base_changes[1],
        base_changes[2],
        base_changes[3],
        base_changes[4],
        base_changes[5],
        base_changes[6],
        base_changes[7],
        base_changes[8],
    )  # i.e. {'TTT': ['A', 'C', 'G'... etc]}

Possible_Amino_Acid_Changes = {}
for k, v in Possible_Codon_Changes.items():
    codon = k
    codon_changes = v
    amino_acid_changes = []
    for i in codon_changes:
        amino_acid_change = Genetic_Code[i]
        amino_acid_changes.append(amino_acid_change)
    Possible_Amino_Acid_Changes[k] = (
        amino_acid_changes[0],
        amino_acid_changes[1],
        amino_acid_changes[2],
        amino_acid_changes[3],
        amino_acid_changes[4],
        amino_acid_changes[5],
        amino_acid_changes[6],
        amino_acid_changes[7],
        amino_acid_changes[8],
    )  # i.e. {'TTT': ['Ile', 'Leu', 'Val'...etc]}

Possible_Amino_Acid_Changes_letters = {}
for k, v in Possible_Amino_Acid_Changes.items():
    codon = k
    codon_changes = v
    amino_acid_changes_letters = []
    for i in codon_changes:
        amino_acid_change_letter = Amino_acid_abbreviations[i]
        amino_acid_changes_letters.append(amino_acid_change_letter)
    Possible_Amino_Acid_Changes_letters[k] = (
        amino_acid_changes_letters[0],
        amino_acid_changes_letters[1],
        amino_acid_changes_letters[2],
        amino_acid_changes_letters[3],
        amino_acid_changes_letters[4],
        amino_acid_changes_letters[5],
        amino_acid_changes_letters[6],
        amino_acid_changes_letters[7],
        amino_acid_changes_letters[8],
    )  # i.e. {'TTT': ['I', 'L', 'V'...etc]}

# Import mutation rates
# csv file containing trinucleotide-context specific mutation rates (calculated from Lee-Six et al 2018 data)
script_dir = os.path.dirname(os.path.realpath(__file__))
filename = os.path.join(script_dir, "Trinucleotide_context_mutation_rates.csv")

with open(filename, "r") as csvfile:
    # csv.reader returns a reader object which will iterate over lines in the csvfile
    read_reader = csv.reader(csvfile)
    row_count = 0
    Mutation_rates = {}

    for row in read_reader:
        if row_count > 0:
            #             print(row)
            site = row[0]
            comp_site = row[2]
            if site != "":
                mutation_rate = float(row[7])
                Mutation_rates[site] = mutation_rate
                Mutation_rates[comp_site] = mutation_rate
        row_count = row_count + 1


# in format (882, 'R882H')
def mutation_rates_list_variants_gene(variant_tuple, gene_exons, gene_exons_extra_bases, intron_triplets):
    # The mutation rate at a base is calculated based on its trinucleotide context.  To obtain the trinucleotide context, the base either side of the mutated base needs to be ascertained.
    # The coding sequence (+ an extra base at each end) of each gene is used when determining the trinucleotide context of the base change.  The problem with using the coding sequence is that
    # if a base change falls at the end of an exon, its neighbouring base will actuall be an intronic base, which may be different to the 'neighbouring' base in the coding sequence.
    # To account for this, there is an "intron_exon_check" function that is called at the intron/exon boundary if the neighbouring base in the intron is different to the neighbouring base in the coding sequence.

    # Split the sequence into codons
    gene_codons = [gene_exons[i : i + 3] for i in range(0, len(gene_exons), 3)]

    # Create a list of the codons they called variants at:
    gene_codons_called = []
    gene_variant_DICT = {}

    codon0 = variant_tuple[0] - 1
    codon1 = variant_tuple[0]
    variant = variant_tuple[1]
    gene_codons_called.append(codon0)  # List of codon numbers that they called variants at
    gene_variant_DICT[row_count] = (
        codon1,
        variant,
    )  # Dictionary:  K = codon number, V = variant at that codon (e.g. 76: p.V76L)

    # Looking at what possible amino acid changes are possible at those codons
    gene_codons_DICT = {}
    calls_gene = {}

    for i in gene_codons_called:  # Iterates over the list of codon numbers
        gene_codons_DICT[i + 1] = gene_codons[
            i
        ]  # Creates of dicitonary of key = codon number, value = bases (sequence) at that codon
        # e.g. 12: CGT, 882: GCA

    for k, v in gene_codons_DICT.items():  # iterating over the list of codon numbers (and codon sequences)
        if (
            v in Possible_Amino_Acid_Changes_letters.keys()
        ):  # assign the list of all the possible amino acid changes to each codon they called
            changes = Possible_Amino_Acid_Changes_letters[v]
            #         print(Possible_Amino_Acid_Changes.values())
            calls_gene[k] = (v, changes)
        # e.g. 12: (CGT, 'K', '*', 'T'etc...)

    # Changing the original amino acid to a Word (e.g. Phe) rather than a codon
    calls_gene_extra = {}

    for k, v in calls_gene.items():
        if v[0] in Genetic_Code:
            Amino_acid_change = Genetic_Code.get(v[0])
            calls_gene_extra[k] = (Amino_acid_change, v[1])
            # e.g. 12: (Phe, 'K', '*', 'T' etc..)

    # Changing the original amino acid to a LETTER (rather than the codon)
    calls_gene_extra2 = {}

    for k, v in calls_gene_extra.items():
        if v[0] in Amino_acid_abbreviations:
            Amino_acid_abbrev = Amino_acid_abbreviations.get(v[0])
            calls_gene_extra2[k] = (Amino_acid_abbrev, v[1])
            # e.g. 12: ((F), ('K', '*', 'T' etc..))

    # Create a dictionary containing info on which position in the codon is affected to create the variant
    gene_position_in_list_dict = {}
    for (
        k,
        v,
    ) in (
        gene_variant_DICT.items()
    ):  # Iterates over the dictionary of (key = row_count, value = aa_pos, variant (e.g. 3: (76, p.V76L))
        observed_change = (v[1])[-1]  # Selects the amino acid change (e.g. the L from p.V76L)
        possible_changes = calls_gene[v[0]][
            1
        ]  # possible changes are the second value in the dictionary of (k = codon number, value = (original amino acid, possible changes))
        possible_changes_info = calls_gene[v[0]]
        codon_seq = calls_gene[v[0]][0]
        possible_changes_string = "".join(possible_changes)
        if (
            v[1][0] == calls_gene_extra2[v[0]][0]
        ):  # v[1][2] = the original amino acid letter (i.e if this aa is in this transcript)
            position_in_list = -1
            positions_in_list = []
            while True:
                position_in_list = possible_changes_string.find(
                    observed_change, position_in_list + 1
                )  # counts position in possible changes (if 0-2 = 1st base of codon, if 3-5 = 2nd base, if 6-8 = 3rd base)
                if position_in_list == -1:
                    break
                positions_in_list.append(position_in_list)
                gene_position_in_list_dict[k] = (positions_in_list, possible_changes_info, codon_seq, v[0], v[1])

    aa_positions_variants_called = []
    for k, v in gene_position_in_list_dict.items():
        aa_positions_variants_called.append(v[3])
    #     print('positions at which nonsense variants called = ', aa_positions_variants_called)

    # Create a dictionary that contains the position in the list of possible AA where the variant is and what it changes to
    gene_bases = {}
    for k, v in gene_position_in_list_dict.items():
        codon_seq = v[2]
        base_changes = []
        if codon_seq in Possible_Base_Changes.keys():
            list_of_possible_bases = Possible_Base_Changes.get(codon_seq)
        for i in v[0]:
            base_change_position = i
            base_change = list_of_possible_bases[base_change_position]
            base_changes.append(base_change)
        gene_bases[k] = (v, base_changes)

    # Create dictionary which contains original codon seq and details on which base is changed
    gene_bases2 = {}

    for k, v in gene_bases.items():
        base_change_positions = v[0][0]
        codon_positions = []
        positions_in_gene = []
        for i in base_change_positions:
            base_change_position = i
            if base_change_position <= 2:
                codon_pos = 1
            if 2 < base_change_position <= 5:
                codon_pos = 2
            if 5 < base_change_position <= 8:
                codon_pos = 3
            codon_positions.append(codon_pos)
        gene_bases2[k] = (v, codon_positions)

    # Create dictionary which contains original codon seq, details on which base is changed and it's position in the coding region
    gene_bases3 = {}

    for k, v in gene_bases2.items():
        codon_positions = v[1]
        base_changes = v[0][1]
        original_bases = []
        positions_in_gene = []
        codon_number = v[0][0][3]
        variant = v[0][0][4]
        for i in codon_positions:
            codon_pos = i
            position_in_gene = (
                (int(v[0][0][3]) - 1) * 3
            ) + codon_pos  # codon number (v[0][0][3]) minus 1 (e.g. codon 307 -1 = 306, then times by 3 takes you to beginning of codon 307 (then add the position))
            positions_in_gene.append(position_in_gene)
            original_base = v[0][0][2][codon_pos - 1]
            original_bases.append(original_base)
        gene_bases3[k] = (v[0][0][1], positions_in_gene, original_bases, base_changes, codon_number, variant)

    # Create dictionary which contains original codon seq, details on which base is changed and its neighbouring bases
    gene_bases4 = {}

    for k, v in gene_bases3.items():
        position_in_gene = v[1]
        sites_trio_nonsense = []

        for i in position_in_gene:
            if i in intron_triplets.keys():
                adapted_trio = intron_triplets.get(i, None)
                sites_trio_nonsense.append([adapted_trio])
                gene_bases4[k] = (v, sites_trio_nonsense)

            else:
                for i in position_in_gene:
                    gene_pos = i
                    site_trio_nonsense = [
                        gene_exons_extra_bases[gene_pos : gene_pos + 3] for gene_pos in range(gene_pos - 1, gene_pos, 3)
                    ]
                    sites_trio_nonsense.append(site_trio_nonsense)
                    gene_bases4[k] = (v, sites_trio_nonsense)

    # Create a dictionary containing pos in coding region, original base, new base, base affected with its neighbours
    # e.g. 938, G, A, TGG = base 938 on coding region, original base = G, changes to base A, base with neighbours either side = TGG
    gene_sites5 = {}

    a = 0
    for k, v in gene_bases4.items():
        codon_n = v[0][1]  # base number in codon (out of 3 x 912)
        codon_n_912 = v[0][4]
        original_base = v[0][2]
        new_base = v[0][3]
        site_trio = v[1]
        info = []
        variant = v[0][5]

        change = 0
        for i in codon_n:
            codon_number = i
            orig_base = original_base[change]
            n_base = new_base[change]
            site_tri = site_trio[change]
            info.append([codon_number, orig_base, n_base, site_tri, codon_n_912, variant])
            change = change + 1
        gene_sites5[a] = info
        a = a + 1
    #     print(gene_sites5)

    # Write in a way that can iterate over mutation rates (i.e. put in to form C[C>A]G etc)
    gene_sites6 = {}

    separating_positions = {}  # separate out all the possible changes to just 1 effective item in the list
    a = 0
    for k, v in gene_sites5.items():
        for i in v:  # because some contain more than 1 change from that codon:
            separating_positions[a] = i
            a = a + 1

    for k, v in separating_positions.items():
        original_base = v[1]
        new_base = v[2]
        site_trio = v[3]
        codon_n_912 = v[4]
        variant = v[5]
        sequence = (str(original_base), str(new_base))
        changes = ">".join(sequence)
        left_base = site_trio[0][0]
        right_base = site_trio[0][2]
        sequence2 = (changes, right_base)
        changes2 = "]".join(sequence2)
        sequence3 = (left_base, changes2)
        changes3 = "[".join(sequence3)
        gene_sites6[k] = (v[0], changes3, codon_n_912, variant)  # site_trio = the affected base and its neighbours
    #         print(gene_sites6)

    # Assign the mutation rates to each site...
    gene_mutation_rates = {}
    for k, v in gene_sites6.items():
        base_change_info = v[1]
        if base_change_info in Mutation_rates.keys():
            mutation_rate = Mutation_rates.get(base_change_info)
            gene_mutation_rates[k] = (v, mutation_rate)

    gene_mutation_rate_sum = []
    for k, v in gene_mutation_rates.items():
        mutation_rate_value = v[1]
        variant = v[0][3]
        gene_mutation_rate_sum.append(mutation_rate_value)

    gene_Summed_mutation_rate_added = np.sum(gene_mutation_rate_sum)

    return gene_Summed_mutation_rate_added


geneConfig = {"DNMT3A": {}}
# gene transcript:
# Downloaded from http://grch37.ensembl.org
geneConfig["DNMT3A"][
    "gene_exons"
] = "ATGCCCGCCATGCCCTCCAGCGGCCCCGGGGACACCAGCAGCTCTGCTGCGGAGCGGGAGGAGGACCGAAAGGACGGAGAGGAGCAGGAGGAGCCGCGTGGCAAGGAGGAGCGCCAAGAGCCCAGCACCACGGCACGGAAGGTGGGGCGGCCTGGGAGGAAGCGCAAGCACCCCCCGGTGGAAAGCGGTGACACGCCAAAGGACCCTGCGGTGATCTCCAAGTCCCCATCCATGGCCCAGGACTCAGGCGCCTCAGAGCTATTACCCAATGGGGACTTGGAGAAGCGGAGTGAGCCCCAGCCAGAGGAGGGGAGCCCTGCTGGGGGGCAGAAGGGCGGGGCCCCAGCAGAGGGAGAGGGTGCAGCTGAGACCCTGCCTGAAGCCTCAAGAGCAGTGGAAAATGGCTGCTGCACCCCCAAGGAGGGCCGAGGAGCCCCTGCAGAAGCGGGCAAAGAACAGAAGGAGACCAACATCGAATCCATGAAAATGGAGGGCTCCCGGGGCCGGCTGCGGGGTGGCTTGGGCTGGGAGTCCAGCCTCCGTCAGCGGCCCATGCCGAGGCTCACCTTCCAGGCGGGGGACCCCTACTACATCAGCAAGCGCAAGCGGGACGAGTGGCTGGCACGCTGGAAAAGGGAGGCTGAGAAGAAAGCCAAGGTCATTGCAGGAATGAATGCTGTGGAAGAAAACCAGGGGCCCGGGGAGTCTCAGAAGGTGGAGGAGGCCAGCCCTCCTGCTGTGCAGCAGCCCACTGACCCCGCATCCCCCACTGTGGCTACCACGCCTGAGCCCGTGGGGTCCGATGCTGGGGACAAGAATGCCACCAAAGCAGGCGATGACGAGCCAGAGTACGAGGACGGCCGGGGCTTTGGCATTGGGGAGCTGGTGTGGGGGAAACTGCGGGGCTTCTCCTGGTGGCCAGGCCGCATTGTGTCTTGGTGGATGACGGGCCGGAGCCGAGCAGCTGAAGGCACCCGCTGGGTCATGTGGTTCGGAGACGGCAAATTCTCAGTGGTGTGTGTTGAGAAGCTGATGCCGCTGAGCTCGTTTTGCAGTGCGTTCCACCAGGCCACGTACAACAAGCAGCCCATGTACCGCAAAGCCATCTACGAGGTCCTGCAGGTGGCCAGCAGCCGCGCGGGGAAGCTGTTCCCGGTGTGCCACGACAGCGATGAGAGTGACACTGCCAAGGCCGTGGAGGTGCAGAACAAGCCCATGATTGAATGGGCCCTGGGGGGCTTCCAGCCTTCTGGCCCTAAGGGCCTGGAGCCACCAGAAGAAGAGAAGAATCCCTACAAAGAAGTGTACACGGACATGTGGGTGGAACCTGAGGCAGCTGCCTACGCACCACCTCCACCAGCCAAAAAGCCCCGGAAGAGCACAGCGGAGAAGCCCAAGGTCAAGGAGATTATTGATGAGCGCACAAGAGAGCGGCTGGTGTACGAGGTGCGGCAGAAGTGCCGGAACATTGAGGACATCTGCATCTCCTGTGGGAGCCTCAATGTTACCCTGGAACACCCCCTCTTCGTTGGAGGAATGTGCCAAAACTGCAAGAACTGCTTTCTGGAGTGTGCGTACCAGTACGACGACGACGGCTACCAGTCCTACTGCACCATCTGCTGTGGGGGCCGTGAGGTGCTCATGTGCGGAAACAACAACTGCTGCAGGTGCTTTTGCGTGGAGTGTGTGGACCTCTTGGTGGGGCCGGGGGCTGCCCAGGCAGCCATTAAGGAAGACCCCTGGAACTGCTACATGTGCGGGCACAAGGGTACCTACGGGCTGCTGCGGCGGCGAGAGGACTGGCCCTCCCGGCTCCAGATGTTCTTCGCTAATAACCACGACCAGGAATTTGACCCTCCAAAGGTTTACCCACCTGTCCCAGCTGAGAAGAGGAAGCCCATCCGGGTGCTGTCTCTCTTTGATGGAATCGCTACAGGGCTCCTGGTGCTGAAGGACTTGGGCATTCAGGTGGACCGCTACATTGCCTCGGAGGTGTGTGAGGACTCCATCACGGTGGGCATGGTGCGGCACCAGGGGAAGATCATGTACGTCGGGGACGTCCGCAGCGTCACACAGAAGCATATCCAGGAGTGGGGCCCATTCGATCTGGTGATTGGGGGCAGTCCCTGCAATGACCTCTCCATCGTCAACCCTGCTCGCAAGGGCCTCTACGAGGGCACTGGCCGGCTCTTCTTTGAGTTCTACCGCCTCCTGCATGATGCGCGGCCCAAGGAGGGAGATGATCGCCCCTTCTTCTGGCTCTTTGAGAATGTGGTGGCCATGGGCGTTAGTGACAAGAGGGACATCTCGCGATTTCTCGAGTCCAACCCTGTGATGATTGATGCCAAAGAAGTGTCAGCTGCACACAGGGCCCGCTACTTCTGGGGTAACCTTCCCGGTATGAACAGGCCGTTGGCATCCACTGTGAATGATAAGCTGGAGCTGCAGGAGTGTCTGGAGCATGGCAGGATAGCCAAGTTCAGCAAAGTGAGGACCATTACTACGAGGTCAAACTCCATAAAGCAGGGCAAAGACCAGCATTTTCCTGTCTTCATGAATGAGAAAGAGGACATCTTATGGTGCACTGAAATGGAAAGGGTATTTGGTTTCCCAGTCCACTATACTGACGTCTCCAACATGAGCCGCTTGGCGAGGCAGAGACTGCTGGGCCGGTCATGGAGCGTGCCAGTCATCCGCCACCTCTTCGCTCCGCTGAAGGAGTATTTTGCGTGTGTGTAA"
geneConfig["DNMT3A"]["gene_exons_extra_bases"] = f"G{geneConfig['DNMT3A']['gene_exons']}G"
geneConfig["DNMT3A"]["intron_exon_boundary_caution_base_gene"] = [
    1279,
    1429,
    1474,
    1554,
    1852,
    2083,
    2082,
    2173,
    2322,
    2478,
]
gene_exons_extra_bases = geneConfig["DNMT3A"]["gene_exons_extra_bases"]
geneConfig["DNMT3A"]["intron_triplets"] = {
    1279: gene_exons_extra_bases[1278] + gene_exons_extra_bases[1279] + "G",
    1429: gene_exons_extra_bases[1428] + gene_exons_extra_bases[1429] + "G",
    1474: gene_exons_extra_bases[1473] + gene_exons_extra_bases[1474] + "G",
    1554: gene_exons_extra_bases[1553] + gene_exons_extra_bases[1554] + "G",
    1852: "G" + gene_exons_extra_bases[1852] + gene_exons_extra_bases[1853],
    2083: "G" + gene_exons_extra_bases[2083] + gene_exons_extra_bases[2084],
    2082: gene_exons_extra_bases[2081] + gene_exons_extra_bases[2082] + "G",
    2173: gene_exons_extra_bases[2172] + gene_exons_extra_bases[2173] + "G",
    2322: gene_exons_extra_bases[2321] + gene_exons_extra_bases[2322] + "G",
    2478: gene_exons_extra_bases[2477] + gene_exons_extra_bases[2478] + "G",
}


if __name__ == "__main__":
    # gene transcript:
    # Downloaded from http://grch37.ensembl.org
    gene_exons = "ATGCCCGCCATGCCCTCCAGCGGCCCCGGGGACACCAGCAGCTCTGCTGCGGAGCGGGAGGAGGACCGAAAGGACGGAGAGGAGCAGGAGGAGCCGCGTGGCAAGGAGGAGCGCCAAGAGCCCAGCACCACGGCACGGAAGGTGGGGCGGCCTGGGAGGAAGCGCAAGCACCCCCCGGTGGAAAGCGGTGACACGCCAAAGGACCCTGCGGTGATCTCCAAGTCCCCATCCATGGCCCAGGACTCAGGCGCCTCAGAGCTATTACCCAATGGGGACTTGGAGAAGCGGAGTGAGCCCCAGCCAGAGGAGGGGAGCCCTGCTGGGGGGCAGAAGGGCGGGGCCCCAGCAGAGGGAGAGGGTGCAGCTGAGACCCTGCCTGAAGCCTCAAGAGCAGTGGAAAATGGCTGCTGCACCCCCAAGGAGGGCCGAGGAGCCCCTGCAGAAGCGGGCAAAGAACAGAAGGAGACCAACATCGAATCCATGAAAATGGAGGGCTCCCGGGGCCGGCTGCGGGGTGGCTTGGGCTGGGAGTCCAGCCTCCGTCAGCGGCCCATGCCGAGGCTCACCTTCCAGGCGGGGGACCCCTACTACATCAGCAAGCGCAAGCGGGACGAGTGGCTGGCACGCTGGAAAAGGGAGGCTGAGAAGAAAGCCAAGGTCATTGCAGGAATGAATGCTGTGGAAGAAAACCAGGGGCCCGGGGAGTCTCAGAAGGTGGAGGAGGCCAGCCCTCCTGCTGTGCAGCAGCCCACTGACCCCGCATCCCCCACTGTGGCTACCACGCCTGAGCCCGTGGGGTCCGATGCTGGGGACAAGAATGCCACCAAAGCAGGCGATGACGAGCCAGAGTACGAGGACGGCCGGGGCTTTGGCATTGGGGAGCTGGTGTGGGGGAAACTGCGGGGCTTCTCCTGGTGGCCAGGCCGCATTGTGTCTTGGTGGATGACGGGCCGGAGCCGAGCAGCTGAAGGCACCCGCTGGGTCATGTGGTTCGGAGACGGCAAATTCTCAGTGGTGTGTGTTGAGAAGCTGATGCCGCTGAGCTCGTTTTGCAGTGCGTTCCACCAGGCCACGTACAACAAGCAGCCCATGTACCGCAAAGCCATCTACGAGGTCCTGCAGGTGGCCAGCAGCCGCGCGGGGAAGCTGTTCCCGGTGTGCCACGACAGCGATGAGAGTGACACTGCCAAGGCCGTGGAGGTGCAGAACAAGCCCATGATTGAATGGGCCCTGGGGGGCTTCCAGCCTTCTGGCCCTAAGGGCCTGGAGCCACCAGAAGAAGAGAAGAATCCCTACAAAGAAGTGTACACGGACATGTGGGTGGAACCTGAGGCAGCTGCCTACGCACCACCTCCACCAGCCAAAAAGCCCCGGAAGAGCACAGCGGAGAAGCCCAAGGTCAAGGAGATTATTGATGAGCGCACAAGAGAGCGGCTGGTGTACGAGGTGCGGCAGAAGTGCCGGAACATTGAGGACATCTGCATCTCCTGTGGGAGCCTCAATGTTACCCTGGAACACCCCCTCTTCGTTGGAGGAATGTGCCAAAACTGCAAGAACTGCTTTCTGGAGTGTGCGTACCAGTACGACGACGACGGCTACCAGTCCTACTGCACCATCTGCTGTGGGGGCCGTGAGGTGCTCATGTGCGGAAACAACAACTGCTGCAGGTGCTTTTGCGTGGAGTGTGTGGACCTCTTGGTGGGGCCGGGGGCTGCCCAGGCAGCCATTAAGGAAGACCCCTGGAACTGCTACATGTGCGGGCACAAGGGTACCTACGGGCTGCTGCGGCGGCGAGAGGACTGGCCCTCCCGGCTCCAGATGTTCTTCGCTAATAACCACGACCAGGAATTTGACCCTCCAAAGGTTTACCCACCTGTCCCAGCTGAGAAGAGGAAGCCCATCCGGGTGCTGTCTCTCTTTGATGGAATCGCTACAGGGCTCCTGGTGCTGAAGGACTTGGGCATTCAGGTGGACCGCTACATTGCCTCGGAGGTGTGTGAGGACTCCATCACGGTGGGCATGGTGCGGCACCAGGGGAAGATCATGTACGTCGGGGACGTCCGCAGCGTCACACAGAAGCATATCCAGGAGTGGGGCCCATTCGATCTGGTGATTGGGGGCAGTCCCTGCAATGACCTCTCCATCGTCAACCCTGCTCGCAAGGGCCTCTACGAGGGCACTGGCCGGCTCTTCTTTGAGTTCTACCGCCTCCTGCATGATGCGCGGCCCAAGGAGGGAGATGATCGCCCCTTCTTCTGGCTCTTTGAGAATGTGGTGGCCATGGGCGTTAGTGACAAGAGGGACATCTCGCGATTTCTCGAGTCCAACCCTGTGATGATTGATGCCAAAGAAGTGTCAGCTGCACACAGGGCCCGCTACTTCTGGGGTAACCTTCCCGGTATGAACAGGCCGTTGGCATCCACTGTGAATGATAAGCTGGAGCTGCAGGAGTGTCTGGAGCATGGCAGGATAGCCAAGTTCAGCAAAGTGAGGACCATTACTACGAGGTCAAACTCCATAAAGCAGGGCAAAGACCAGCATTTTCCTGTCTTCATGAATGAGAAAGAGGACATCTTATGGTGCACTGAAATGGAAAGGGTATTTGGTTTCCCAGTCCACTATACTGACGTCTCCAACATGAGCCGCTTGGCGAGGCAGAGACTGCTGGGCCGGTCATGGAGCGTGCCAGTCATCCGCCACCTCTTCGCTCCGCTGAAGGAGTATTTTGCGTGTGTGTAA"

    # gene transcript:
    # Downloaded from http://grch37.ensembl.org
    gene_exons_extra_bases = f"G{gene_exons}G"

    intron_exon_boundary_caution_base_gene = [1279, 1429, 1474, 1554, 1852, 2083, 2082, 2173, 2322, 2478]

    intron_triplets = {
        1279: gene_exons_extra_bases[1278] + gene_exons_extra_bases[1279] + "G",
        1429: gene_exons_extra_bases[1428] + gene_exons_extra_bases[1429] + "G",
        1474: gene_exons_extra_bases[1473] + gene_exons_extra_bases[1474] + "G",
        1554: gene_exons_extra_bases[1553] + gene_exons_extra_bases[1554] + "G",
        1852: "G" + gene_exons_extra_bases[1852] + gene_exons_extra_bases[1853],
        2083: "G" + gene_exons_extra_bases[2083] + gene_exons_extra_bases[2084],
        2082: gene_exons_extra_bases[2081] + gene_exons_extra_bases[2082] + "G",
        2173: gene_exons_extra_bases[2172] + gene_exons_extra_bases[2173] + "G",
        2322: gene_exons_extra_bases[2321] + gene_exons_extra_bases[2322] + "G",
        2478: gene_exons_extra_bases[2477] + gene_exons_extra_bases[2478] + "G",
    }

    print(
        mutation_rates_list_variants_gene(
            [882, "R882C"], gene_exons, gene_exons_extra_bases, intron_triplets
        )  # 1.41512e-08
    )
    print(
        mutation_rates_list_variants_gene(
            [882, "R882H"], gene_exons, gene_exons_extra_bases, intron_triplets
        )  # 1.88229e-08
    )

# %%
