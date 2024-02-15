import subprocess
import sys

weed_cmd_rep = "ska weed -o no_ambig.skf --min-freq 1 --filter no-ambig-or-const test_grid.skf"
weed_cmd_norep = "ska weed -o ambig.skf --min-freq 1 --filter no-const test_grid.skf"

def compare_hits(repeats: bool):
    found = set()
    expected = set()
    if repeats:
        nk_cmd = "ska nk --full-info no_ambig.skf | awk 'NR > 9 {print $0}' | sed '$ d' > test1.txt"
        subprocess.run(nk_cmd, shell=True, stderr=subprocess.DEVNULL)
        f = open("test1.txt", 'r')
        for line in f:
            found.add(line.rstrip())
        f.close()

        f = open("split_kmers.txt", 'r')
        for line in f:
            expected.add(line.rstrip())
        f.close()
    else:
        nk_cmd = "ska nk --full-info ambig.skf | awk 'NR > 9 {print $0}' | sed '$ d' > test1.txt"
        subprocess.run(nk_cmd, shell=True, stderr=subprocess.DEVNULL)
        f = open("test1.txt", 'r')
        for line in f:
            found.add(str(line.rstrip().split("\t")[0:2]))
        f.close()

        f = open("split_kmers.txt", 'r')
        for line in f:
            expected.add(str(line.rstrip().split("\t")[0:2]))
        f.close()
    power = len(found.intersection(expected)) / len(expected)
    return power

print(f"k\tpi\tindel_rate\tRepeats\tPower")
for k in [17, 31, 63]:
    for dist in [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1]:
        for indel_rate in [0.05, 0.1, 0.2]:
            for repeat in range(10):
                sim_cmd = f"ska_sim -k {k} -e {dist} -i {indel_rate} > seq.fa"
                subprocess.run(sim_cmd, shell=True)
                ska_cmd = f"ska build --single-strand -k {k} -o test_grid pneumo.fa seq.fa"
                subprocess.run(ska_cmd, shell=True, stderr=subprocess.DEVNULL)

                subprocess.run(weed_cmd_rep, shell=True, stderr=subprocess.DEVNULL)
                overlap = compare_hits(True)
                print(f"{k}\t{dist}\t{indel_rate}\t{repeat}\tNo ambiguity\t{overlap}")

                subprocess.run(weed_cmd_norep, shell=True, stderr=subprocess.DEVNULL)
                overlap = compare_hits(False)
                print(f"{k}\t{dist}\t{indel_rate}\t{repeat}\tAllow ambiguity\t{overlap}")

sys.exit(0)