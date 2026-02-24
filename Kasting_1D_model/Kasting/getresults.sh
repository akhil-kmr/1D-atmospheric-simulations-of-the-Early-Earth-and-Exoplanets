#!/bin/bash
#!/bin/awk
#set echo
#set verbose

#############
# WROTE BY PENG LIU AT JUL/2019 IN PENN STATE
# THIS SCRIPT CAN HELP YOU FIND DATA YOU WANT FROM OUTPUTS
#
#############

OXYR=7 # YOU GOT 7 O2 MIXING RATIO
OXYPAL=(0.0 1.0 0.3 0.1 0.03 0.01 0.003 0.001) # O2 MIXING RATIO RANGE [1PAL-0.001PAL]
CO2R=6 # YOU GOT 6 CO2 MIXING RATIO
CO2PAL=(0.0 1.0 3.0 10.0 30.0 100.0 300.0)    # CO2 MIXING RATIO RANGE [1PAL-300PAL]

#DATA LOGIC
#     O2 LEVELS FOR 6 CO2 LEVELS
#NEXT O2 LEVELS FOR 6 CO2 LEVELS

###############################
#GET GPP FOR BOTH O2 AND CO2
##############################
echo "GPP  O2  CO2" >>data_plot.dat
i=1
j=1
while [ $i -le ${OXYR} ];do
j=1
while [ $j -le ${CO2R} ];do
row1=4
cat ./IO/history/O2_${OXYPAL[i]}PAL_CO2_${CO2PAL[j]}PAL/isotope.dat | awk 'NR=='''$row1''' {print $7,$8}' >>data_plot.dat
  j=$((j+1))
done

  i=$((i+1))
done

###############################
#GET DELTA FOR BOTH O2 AND H2SO4 CO2 AT THE SURFACE
##############################
echo ' '
echo "delta  O2  CO2 H2SO4" >>data_plot.dat
row2=900
i=1
j=1
while [ $i -le ${OXYR} ];do
j=1
while [ $j -le ${CO2R} ];do
cat ./IO/history/O2_${OXYPAL[i]}PAL_CO2_${CO2PAL[j]}PAL/isotope.dat | awk 'NR=='''$row2''' {print $5,$10,$8}' >>data_plot.dat
  j=$((j+1))
done

  i=$((i+1))
done

###############################
#ALL TERMS ABOUT PRODUCTION AND LOSS
##############################
echo ' '
echo "OXYDEPI OXYRANI OXYUPI TPI-TLI(H2Q) FLOW(LCQ) FLOW(LN2Q)" >>data_plot.dat
row3=1108
row4=1109
i=1
j=1
while [ $i -le ${OXYR} ];do
j=1
while [ $j -le ${CO2R} ];do
cat ./IO/history/O2_${OXYPAL[i]}PAL_CO2_${CO2PAL[j]}PAL/isotope.dat | awk 'NR=='''$row3''' {print $3,$6,$9,$12}' >>data_plot.dat
cat ./IO/history/O2_${OXYPAL[i]}PAL_CO2_${CO2PAL[j]}PAL/isotope.dat | awk 'NR=='''$row4''' {print $3,$6}' >>data_plot.dat
  j=$((j+1))
done

  i=$((i+1))
done


###############################
#TP(H2O)-TL(H2O) VS FO2 for PRESENT CO2 LEVEL
##############################
echo ' '
echo "TP(H2O)-TL(H2O) XPALO2 1PALCO2" >>data_plot.dat
j=1
while [ $j -le ${CO2R} ];do

grep "TP-TL(H2O) =" ./IO/history/O2_${OXYPAL[1]}PAL_CO2_${CO2PAL[j]}PAL/outchem.dat >>data_plot.dat

  j=$((j+1))
done





























