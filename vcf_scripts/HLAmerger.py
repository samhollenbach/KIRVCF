import csv
import VCFparser
import copy



def write_outfile(HLA_alles):
    header = copy.copy(header_master)
    if 'ALL' in HLA_alles:
        name = 'trios_HLA_all.csv'
        indeces = list(range(6,len(header)))
    else:
        indeces = [header.index(i) for i in HLA_alles]
        HLA_alles = [l.replace('HLA-', '') for l in HLA_alles]
        name = 'trios_HLA_{}.csv'.format('_'.join(HLA_alles))

    colname = '~'.join(HLA_alles)
    header = header[:6] + [colname]*2
    print(header)
    #rows = copy.copy(rows_master)
    rows_new = []
    for i in range(0, len(rows_master), 2):
        row = rows_master[i]
        row2 = rows_master[i+1]
        if '' in row:
            continue
        lead = row[:6]
        col1 = "~".join([row[i] for i in indeces])
        col2 = "~".join([row2[i] for i in indeces])
        rows_new.append(lead + [col1] + [col2])

    VCFparser.write_output_file(name, header, rows_new, delim=",", output_dir='../results/trios/')


def oneline_merge(infile):
    with open(infile, 'r') as r:
        reader = csv.reader(r, delimiter=",")
        header_master = next(reader)
        header = header_master[:6]
        data = []
        for h in header_master[6:]:
            header.append(h)
            header.append(h)

        for row in reader:
            row2 = next(reader)
            row_new = row[:6]
            for i, e in enumerate(row[6:], 6):
                row_new.append(row[i])
                row_new.append(row2[i])
            data.append(row_new)
        VCFparser.write_output_file("oneline_{}.txt".format(infile.split("/")[-1].split(".")[0]), header, data, delim="\t",output_dir="../results/")






with open("../results/INDIGO_Haplotype_working_2.csv", 'r') as r:
    reader = csv.reader(r, delimiter=",")
    header_master = next(reader)
    rows_master = [r for r in reader]
    print("MASTER HEADER : ", header_master)
    # write_outfile(['ALL']) #ALL
    # write_outfile(['HLA-DRB1', 'HLA-DQA1'])
    # write_outfile(['HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1'])
    # write_outfile(['HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1'])
    # write_outfile(['HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPB1'])
    # write_outfile(['HLA-B', 'HLA-DRB1'])
    # write_outfile(['HLA-C', 'HLA-B', 'HLA-DRB1'])
    # write_outfile(['HLA-A', 'HLA-C', 'HLA-B', 'HLA-DRB1'])
    # write_outfile(['HLA-DRB1'])
    # write_outfile(['HLA-DQA1'])
    # write_outfile(['HLA-DQB1'])
    # write_outfile(['HLA-DPB1'])
    # write_outfile(['HLA-DPA1'])
    # write_outfile(['HLA-B'])
    # write_outfile(['HLA-C'])
    # write_outfile(['HLA-A'])
    oneline_merge("../results/INDIGO_Haplotype_working.csv")
    #write_outfile(['HLA-DRB1', 'HLA-DRB345'])


