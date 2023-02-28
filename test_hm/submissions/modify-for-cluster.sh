#!/bin/bash

# Copy contents of real_o3.sh into atlas_real_o3.sh
cp ../real_o3.sh atlas_real_o3.sh

# Uncomment the WD line and comment out the other
sed -i.bak 's|#WD=/work/cameron.mills/projects/pycbc/test_hm|WD=/work/cameron.mills/projects/pycbc/test_hm|g' atlas_real_o3.sh
sed -i.bak 's|WD=/Users/camill/projects/pycbc/test_hm|#WD=/Users/camill/projects/pycbc/test_hm|g' atlas_real_o3.sh

# Delete all lines starting with FRAME_COMMAND and add a new one
sed -i.bak '/^FRAME_COMMAND/d' atlas_real_o3.sh
sed -i.bak 's|#FRAME_COMMAND="--frame-type H1:${h1frame} L1:${l1frame} V1:${v1frame}"|FRAME_COMMAND="--frame-type H1:${h1frame} L1:${l1frame} V1:${v1frame}"|g' atlas_real_o3.sh
#head -n 2 atlas_real_o3.sh > temp.txt
#echo 'FRAME_COMMAND="--frame-type H1:${h1frame} L1:${l1frame} V1:${v1frame}"' >> temp.txt
#tail -n +3 atlas_real_o3.sh >> temp.txt
#mv temp.txt atlas_real_o3.sh


# Modify processing scheme
sed -i.bak 's/--processing-scheme [[:alnum:]]*/--processing-scheme mkl/g' atlas_real_o3.sh

