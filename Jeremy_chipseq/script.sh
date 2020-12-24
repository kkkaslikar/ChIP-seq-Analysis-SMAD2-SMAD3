# 
# macs2 predictd -i ./bamfiles/LAP-untreated_Smad3.bam |tee -a predictd.txt
# macs2 predictd -i ./bamfiles/MDA-untreated_Smad3.bam |tee -a predictd.txt
# macs2 predictd -i ./bamfiles/FCC6CLFACXX-wHAPSI016408-113_L5_1.bam |tee -a predictd.txt

# macs2 callpeak -B -t bamfiles/MDA-untreated_Smad3.bam -c bamfiles/input.bam -n untreated/untreated-native-smad3 --nomodel --extsize 217

# macs2 callpeak -B -t bamfiles/FCC6CLFACXX-wHAPSI016408-113_L5_1.bam -c bamfiles/input.bam -n untreated/untreated-lap-smad2 --nomodel --extsize 217

# macs2 callpeak -B -t bamfiles/LAP-untreated_Smad3.bam -c bamfiles/input.bam -n untreated/untreated-lap-smad3 --nomodel --extsize 217

macs2 bdgdiff --t1 ./untreated/untreated-lap-smad2/untreated-lap-smad2_treat_pileup.bdg --c1 ./untreated/untreated-lap-smad2/untreated-lap-smad2_control_lambda.bdg --t2 ./untreated/untreated-lap-smad3/untreated-lap-smad3_treat_pileup.bdg --c2 ./untreated/untreated-lap-smad3/untreated-lap-smad3_control_lambda.bdg --d1 31588061 --d2 31588061 --o-prefix untreated-diff_lap-smad2_vs_lap-smad3

macs2 bdgdiff --t1 ./untreated/untreated-lap-smad2/untreated-lap-smad2_treat_pileup.bdg --c1 ./untreated/untreated-lap-smad2/untreated-lap-smad2_control_lambda.bdg --t2 ./untreated/untreated-native-smad3/untreated-native-smad3_treat_pileup.bdg --c2 ./untreated/untreated-native-smad3/untreated-native-smad3_control_lambda.bdg --d1 31588061 --d2 13426831 --o-prefix ./untreated/untreated-diffpeak_results/untreated-diff_lap-smad2_vs_native-smad3
