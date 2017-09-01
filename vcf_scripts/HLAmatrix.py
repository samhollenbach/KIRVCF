import csv
import VCFparser

trans_files = VCFparser.get_files_in_dir(directory="../results/trans")
trio_files = VCFparser.get_files_in_dir(directory="../results/trios")

include_labels = False

for trans_file, trio_file in zip(trans_files, trio_files):
    name = trio_file.split("/")[-1].replace("trios", "matrix")
    matrix = []
    with open(trans_file, 'r') as t:
        treader = csv.reader(t, delimiter=",")
        next(treader)
        haps = [line[0] for line in treader]
        print(haps)
        if include_labels:
            matrix.append(["TRANS/NOT_TRANS"] + haps)
        else:
            hap_h = [[h] for h in haps]
            VCFparser.write_output_file("labels_{}".format(name), col_names=[], data=hap_h, output_dir="../results/matrix/")

        l = len(haps)
        for h in haps:
            if include_labels:
                matrix.append([h] + [0] * l)
            else:
                matrix.append([0] * l)

        #print(matrix)
        with open(trio_file, 'r') as r:

            reader = csv.reader(r, delimiter=",")
            header = next(reader)

            for row in reader:
                row2 = next(reader)
                row3 = next(reader)

                child_i = [row[5], row2[5], row3[5]].index('2')
                rows = [row[-2:], row2[-2:], row3[-2:]]
                row_child = rows.pop(child_i)

                trans_dict = {"N0":"", "T0": "", "N1":"","T1":""}

                row0 = rows[0]
                row1 = rows[1]
                for h in row0:
                    if h not in row_child:
                            trans_dict["N0"] = row0.pop(row0.index(h))

                for h in row1:
                    if h not in row_child:
                            trans_dict["N1"] = row1.pop(row1.index(h))

                #print("CHILD: ", row_child)
                if len(row0) is 2 and len(row1) is 2:

                    if row0[0] == row0[1]:
                        t0 = row0.pop(0)
                        trans_dict["T0"] = t0
                        trans_dict["N0"] = row0.pop(0)
                        row_child.pop(row_child.index(t0))
                        for h in row1:
                            if h in row_child:
                                trans_dict["T1"] = row1.pop(row1.index(h))
                                trans_dict["N1"] = row1.pop(0)
                                row_child = []

                    elif row1[0] == row1[1]:
                        t1 = row1.pop(0)
                        trans_dict["T1"] = t1
                        trans_dict["N1"] = row1.pop(0)
                        row_child.pop(row_child.index(t1))

                        for h in row0:
                            if h in row_child:
                                trans_dict["T0"] = row0.pop(row0.index(h))
                                trans_dict["N0"] = row0.pop(0)
                                row_child = []
                    else:
                        # Just choose first one for
                        t0 = row0.pop(0)
                        trans_dict["T0"] = t0
                        trans_dict["N0"] = row0.pop(0)
                        row_child.pop(row_child.index(t0))

                        for h in row1:
                            if h in row_child:
                                trans_dict["T1"] = row1.pop(row1.index(h))
                                trans_dict["N1"] = row1.pop(0)
                                row_child = []


                elif len(row0) is 1 and len(row1) is 1:
                    trans_dict["T0"] = row0.pop(0)
                    trans_dict["T1"] = row1.pop(0)
                    row_child = []
                elif len(row0) is 1:
                    t0 = row0.pop(0)
                    trans_dict["T0"] = t0
                    row_child.pop(row_child.index(t0))

                    if row1[0] == row1[1]:
                        t1 = row1.pop(0)
                        trans_dict["T1"] = t1
                        trans_dict["N1"] = row1.pop(0)
                        row_child.pop(row_child.index(t1))
                    else:
                        for h in row1:
                            if h in row_child:
                                t1 = row1.pop(row1.index(h))
                                trans_dict["T1"] = t1
                                trans_dict["N1"] = row1.pop(0)
                                row_child.pop(row_child.index(h))


                elif len(row1) is 1:
                    t1 = row1.pop(0)
                    trans_dict["T1"] = t1
                    row_child.pop(row_child.index(t1))

                    if row0[0] == row0[1]:
                        t0 = row0.pop(0)
                        trans_dict["T0"] = t0
                        trans_dict["N0"] = row0.pop(0)
                        row_child.pop(row_child.index(t0))
                    else:
                        for h in row0:
                            if h in row_child:
                                t0 = row0.pop(row0.index(h))
                                trans_dict["T0"] = t0
                                trans_dict["N0"] = row0.pop(0)
                                row_child.pop(row_child.index(h))


                i_t0 = 0
                i_t1 = 0
                i_n0 = 0
                i_n1 = 0
                for lab, hap in trans_dict.items():

                    i = haps.index(hap)
                    print(i)
                    if include_labels:
                        i += 1
                    if lab == "N0":
                        i_n0 = i
                    elif lab == "N1":
                        i_n1 = i
                    elif lab == "T0":
                        i_t0 = i
                    elif lab == "T1":
                        i_t1 = i

                matrix[i_t0][i_n0] += 1
                matrix[i_t1][i_n1] += 1


        VCFparser.write_output_file(name, col_names=matrix[0], data=matrix[1:], delim=",",
                                    output_dir="../results/matrix/")
