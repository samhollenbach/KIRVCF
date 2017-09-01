import csv
import VCFparser


files = VCFparser.get_files_in_dir(directory="../results/trios")

for file in files:
    name = file.split("/")[-1].replace("trios", "transmission")
    print(name)

    with open(file, 'r') as r:

        trans_head = ["Haplotype", "Transmitted", "Not Transmitted"]
        reader = csv.reader(r, delimiter=",")
        header = next(reader)
        print(file)
        trans_data = []

        for row in reader:
            row2 = next(reader)
            row3 = next(reader)


            id = [row[0], row2[0], row3[0]]
            print(id)


            child_i = [row[5], row2[5], row3[5]].index('2')
            rows = [row[-2:], row2[-2:], row3[-2:]]
            row_child = rows.pop(child_i)
            #rows = [r2 for r1 in rows for r2 in r1]

            #trans_dict = {"N0":"", "T0": "", "N1":"","T1":""}
            trans_dict = {}

            #print(rows)
            #print(row_child)

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


            print(trans_dict)

            for lab, hap in trans_dict.items():
                hap_index = 0
                present = False
                for i, t in enumerate(trans_data):
                    if hap == t[0]:
                        present = True
                        hap_index = i
                        break
                if not present:
                    trans_data.append([hap, 0, 0])
                    hap_index = len(trans_data) - 1

                if lab is "T0" or lab is "T1":
                    trans_data[hap_index][1] += 1
                elif lab is "N0" or lab is "N1":
                    trans_data[hap_index][2] += 1

        VCFparser.write_output_file(name, col_names=trans_head, data=trans_data, delim=",",output_dir="../results/trans/")





