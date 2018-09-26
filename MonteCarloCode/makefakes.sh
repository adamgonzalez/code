#!/bin/bash

model=EPIC_mos2_cflux-flux-sep_fakfit.xcm
#background=PN_IZW1_C_bg.pha
rmf=mos2.rsp
#arf=PN_IZW1.arf

for i in {1..1000}
do
echo "${i} iteration"
# xspec <<EOF
xspec > /dev/null 2>&1 <<EOF
chatter 1
data mos2_opt.pha
@${model}
energies extend low 0.1
energies extend high 100.0
data none
fakeit none
${rmf}
${arf}
y

fakespec${i}.fak
262588.0
exit
EOF
ftgrouppha outfile=fakespec${i}_grp.fak infile=fakespec${i}.fak grouptype=opt respfile=${rmf}
done
