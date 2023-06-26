grep -aob 'GGGCAGCATTCCCATCTCCAGCTCCCCCTCCCCCGACTCGATTGTACGTG' Wt1os.ko.fa

awk -F "\t" '{
if(($4 > 24090)) {  # condition check
     $4=$4+50  # print desired output
  }
  else{$4
  }
  ; 
if(($5 > 24090)) {  # condition check
     $5=$5+50  # print desired output
  }
  else{$5
  } ; print}'  OFS="\t" Wt1os.mod.gtf > Wt1os.ko.mod.gtf