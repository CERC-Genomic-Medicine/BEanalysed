#!/bin/bash
# Render the structural_cluster_analysis.py output in ChimeraX.
# Usage:
#   for x in *_clusters.defattr ; do bash chimeraX_clusters.sh [Your_Model.cif] $x ; done
#
# Opens the 3D model, loads the cluster_id attribute file, greys out the
# non-clustered residues, colours each cluster with a distinct categorical
# palette, shows the clustered side chains, and saves a ChimeraX session
# plus a glTF export next to the attribute file.
CIF_FILE="$1"
BASENAME=$(basename "$CIF_FILE" .cif)
Attribute="$2"
OUTPUT_FILE=$(basename "$Attribute" _clusters.defattr)_clusters
# Run ChimeraX in headless mode
ChimeraX --nogui --cmd " \
open $CIF_FILE; \
open $Attribute; \
preset publication 1; \
color lightgray; \
color byattribute cluster_id palette Paired-12 novalue lightgray key true; \
select ::cluster_id>=1; \
show sel atoms; \
style sel stick; \
~select; \
2dlab text 'Cluster';\
set bg white;\
light flat ; \
key labelSide left/top ;\
key pos 0.06,0.1 size 0.3,0.04;\
2dlabels xpos 0.021 ypos 0.021 ;\
save $OUTPUT_FILE format session; \
graphics quality ribbonSides 6 ribbonDivisions 10; \
save $OUTPUT_FILE.glb format gltf; \
view matrix ; \
exit
"
