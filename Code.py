#install pycirclize
import matplotlib.pyplot as plt
from pycirclize import Circos
from pycirclize.parser import Gff

# Load GFF file
gff_file = ("enterobacteria.gff")
gff = Gff(gff_file)

circos = Circos(sectors={gff.name: gff.range_size})
circos.text("Enterobacteria phage\n(NC_000902)", size=15)

sector = circos.sectors[0]
cds_track = sector.add_track((90, 100))
cds_track.axis(fc="#EEEEEE", ec="none")
# Plot forward CDS
cds_track.genomic_features(
    gff.extract_features("CDS", target_strand=1),
    plotstyle="arrow",
    r_lim=(95, 100),
    fc="salmon",
)
# Plot reverse CDS
cds_track.genomic_features(
    gff.extract_features("CDS", target_strand=-1),
    plotstyle="arrow",
    r_lim=(90, 95),
    fc="skyblue",
)
# Extract CDS product labels
pos_list, labels = [], []
for f in gff.extract_features("CDS"):
    start, end = int(str(f.location.end)), int(str(f.location.start))
    pos = (start + end) / 2
    label = f.qualifiers.get("product", [""])[0]
    if label == "" or label.startswith("hypothetical"):
        continue
    if len(label) > 20:
        label = label[:20] + "..."
    pos_list.append(pos)
    labels.append(label)
# Plot CDS product labels on outer position
cds_track.xticks(
    pos_list,
    labels,
    label_orientation="vertical",
    show_bottom_line=True,
    label_size=6,
    line_kws=dict(ec="grey"),
)
# Plot xticks & intervals on inner position
cds_track.xticks_by_interval(
    interval=5000,
    outer=False,
    show_bottom_line=True,
    label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
    label_orientation="vertical",
    line_kws=dict(ec="grey"),
)

fig = circos.plotfig()

plt.show()

