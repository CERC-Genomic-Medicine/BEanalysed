#!/bin/bash
#for x in *_Sig.defattr ; do bash chimeraX.sh [Your_Model.cif] $x  [Dupp_Strategy];done
CIF_FILE="$1"
BASENAME=$(basename "$CIF_FILE" .cif)
Attribute="$2"
Dupp_Strategy="$3"
OUTPUT_FILE=$(basename "$Attribute" _Sig.defattr)
# Run ChimeraX in headless mode
ChimeraX --nogui --cmd " \
open $CIF_FILE; \
open $Attribute; \
preset publication 1; \
color  byattr $Dupp_Strategy palette RdYlBu novalue gray key true; \
2dlab text 'Log Fold Change';\
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
