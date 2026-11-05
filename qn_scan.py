import urllib.request
import matplotlib.pyplot as plt

URL = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome%3AUP000002311"
urllib.request.urlretrieve(URL, "yeast.fasta")


def read_fasta(path):
    proteins = {}
    name = None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line.split("|")[2].split()[0]
            proteins[name] = ""
        else:
            proteins[name] += line
    return proteins


def qn_profile(seq, w=40):
    profile = []
    for i in range(len(seq) - w + 1):
        window = seq[i:i + w]
        frac = (window.count("Q") + window.count("N")) / w
        profile.append(frac)
    return profile


proteins = read_fasta("yeast.fasta")
for name in proteins:
    if "ERF3" in name or "SUP35" in name:
        print("нашёл:", name)

for k in ["ERF3_YEAST", "URE2_YEAST", "RNQ1_YEAST", "SWI1_YEAST", "MOT3_YEAST"]:
    print(k, "есть в файле" if k in proteins else "НЕТ в файле")
print(len(proteins), "белков")

scores = {}
for name, seq in proteins.items():
    if len(seq) >= 40:
        scores[name] = max(qn_profile(seq))

top = sorted(scores.items(), key=lambda x: x[1], reverse=True)

print("\nТоп-20:")
for name, score in top[:20]:
    print(name, round(score, 2))

print("\nИзвестные прионы:")
known = ["ERF3_YEAST", "URE2_YEAST", "RNQ1_YEAST", "SWI1_YEAST", "MOT3_YEAST", "CYC8_YEAST"]
for k in known:
    if k in scores:
        place = [n for n, s in top].index(k) + 1
        print(k, "доля:", round(scores[k], 2), "место:", place, "из", len(top))

plt.hist(scores.values(), bins=40)
plt.xlabel("Максимальная доля Q+N в окне 40 а.о.")
plt.ylabel("Число белков")
plt.title("Q/N-богатость протеома S. cerevisiae")

for k in known:
    if k in scores:
        plt.axvline(scores[k], color="red", linestyle="--", linewidth=1)
        plt.text(scores[k], 500, k.split("_")[0], rotation=90, color="red", fontsize=8)

plt.savefig("hist.png", dpi=150)
plt.show()

p = qn_profile(proteins["SUP35_YEAST"])
plt.plot(p)
plt.xlabel("Позиция в белке")
plt.ylabel("Доля Q+N")
plt.title("Sup35")
plt.savefig("sup35.png", dpi=150)
plt.show()
known = ["ERF3_YEAST", "URE2_YEAST", "RNQ1_YEAST", "SWI1_YEAST", "MOT3_YEAST"]
p = qn_profile(proteins["ERF3_YEAST"])
print("\nМесто Sup35:", [n for n, s in top].index("ERF3_YEAST") + 1, "из", len(top))
print("Белков с долей Q+N > 0.4:", sum(1 for s in scores.values() if s > 0.4))
print("Белков с долей Q+N > 0.5:", sum(1 for s in scores.values() if s > 0.5))
