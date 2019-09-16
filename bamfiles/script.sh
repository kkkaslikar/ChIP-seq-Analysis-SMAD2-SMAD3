# macs2 callpeak -B -t MDA-treated_Smad3.bam -c input.bam -n smad3 --nomodel --extsize 218

# macs2 callpeak -B -t LAP-treated_Smad3.bam -c input.bam -n lap-smad3 --nomodel --extsize 218

# macs2 callpeak -B -t FCC6CLFACXX-wHAPSI016409-133_L5_1.bam -c input.bam -n lap-smad2 --nomodel --extsize 218

# macs2 bdgdiff --t1 lap-smad2/lap-smad2_treat_pileup.bdg --c1 lap-smad2/lap-smad2_control_lambda.bdg --t2 lap-smad3/lap-smad3_treat_pileup.bdg --c2 lap-smad3/lap-smad3_control_lambda.bdg --d1 31588061 --d2 10967473 --o-prefix diff_lap-smad2_vs_lap-smad3

macs2 bdgdiff --t1 lap-smad2/lap-smad2_treat_pileup.bdg --c1 lap-smad2/lap-smad2_control_lambda.bdg --t2 native-smad3/native-smad3_treat_pileup.bdg --c2 native-smad3/native-smad3_control_lambda.bdg --d1 31588061 --d2 19110302 --o-prefix diff_lap-smad2_vs_native-smad3
