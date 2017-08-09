import csv
import os

input_dir = "../results/vars"
file_ending = "_variants.txt"
baseline_country = "Germany_part1"

variant_filenames = [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith(file_ending)]

use_baseline = False


def find_baseline(country):
    for file in variant_filenames:
        if country in file:
            return file

    print("Baseline '{}' could not be found".format(country))



def read_singletons(filename, is_baseline = False):
    with open(filename, "r") as r:

        print("Reading {} for singletons".format(filename))

        reader = csv.reader(r, delimiter="\t")
        country = filename.split("/")[-1].replace(file_ending,"")
        print(country)
        singletons = []
        for row in reader:
            if not row:
                continue
            if "Multi Variant Files" in row[0]:
                break
            if "Position:" not in row[0]:
                continue

            row = {r[0] : r[1] for r in map(lambda x: x.split(": "),row)}
            if is_baseline:
                row["Countries"] = [country]

            singletons.append(row)

        return country, singletons


def create_variants_obj():
    baseline_file = find_baseline(baseline_country)

    _, baseline = read_singletons(baseline_file, True)



    for filename in variant_filenames :
        if baseline_country in filename:
            continue
        baseline_positions = [rb['Position'] for rb in baseline]
        country, singletons = read_singletons(filename)
        for s in singletons:
            if s['Position'] in baseline_positions:
                for d in baseline:
                    if d['Position'] == s['Position']:
                        d['Countries'].append(country)
            elif not use_baseline:
                s["Countries"] = [country]
                baseline.append(s)


    return baseline

vars = create_variants_obj()
print(vars)
vars_counts = {}

for v in vars:
    countries_list = v['Countries']
    shares = len(countries_list) - 1
    if shares not in vars_counts:
        vars_counts[shares] = {}
    for c in countries_list:
        if c not in vars_counts[shares]:
            vars_counts[shares][c] = 0
        vars_counts[shares][c] += 1

print(vars_counts)

wd = list(vars_counts.keys())
ws = list(vars_counts.values())



with open("../results/singleton_vars_all.txt", 'w') as w:
    writer = csv.writer(w, delimiter="\t")
    writer.writerow(["Shares", "Country", "Count"])
    for shares, country_counts in vars_counts.items():
        for country, count in country_counts.items():
            writer.writerow([shares, country, count])








