import subprocess
import sys
from tqdm import tqdm

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
    if len(expected) > 0:
        power = len(found.intersection(expected)) / len(expected)
    else:
        power = 1.0
    return power

def avg_range(overlaps: list):
    overlaps.sort()
    overlaps = overlaps[1:20:]
    avg = sum(overlaps) / len(overlaps)
    top = max(overlaps)
    bottom = min(overlaps)
    return (bottom, avg, top)

def main():
    print(f"k\tpi\tindel_rate\tPower_b\tAverage power\tPower_t")
    with tqdm(total=3*7*3*20) as pbar:
        for k in [17, 31, 63]:
            for dist in [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.03, 0.1]:
                for indel_rate in [0.0, 0.1, 0.2]:
                    overlap_rep = list()
                    overlap_norep = list()
                    for _repeat in range(20):
                        sim_cmd = f"ska_sim -k {k} -e {dist} -i {indel_rate} > seq.fa"
                        subprocess.run(sim_cmd, shell=True)
                        ska_cmd = f"ska build --single-strand -k {k} -o test_grid pneumo.fa seq.fa"
                        subprocess.run(ska_cmd, shell=True, stderr=subprocess.DEVNULL)

                        subprocess.run(weed_cmd_rep, shell=True, stderr=subprocess.DEVNULL)
                        overlap_rep.append(compare_hits(True))
                        subprocess.run(weed_cmd_norep, shell=True, stderr=subprocess.DEVNULL)
                        overlap_norep.append(compare_hits(False))
                        pbar.update(1)

                    stats = avg_range(overlap_rep)
                    print(f"{k}\t{dist}\t{indel_rate}\tNo ambiguity\t{stats[0]}\t{stats[1]}\t{stats[2]}")
                    stats = avg_range(overlap_norep)
                    print(f"{k}\t{dist}\t{indel_rate}\tAllow ambiguity\t{stats[0]}\t{stats[1]}\t{stats[2]}")


    sys.exit(0)

if __name__ == "__main__":
    main()
