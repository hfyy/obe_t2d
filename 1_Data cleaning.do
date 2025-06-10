/* this do file is used to create an analytic file from RAND longitudinal HRS 1992-2020*/
clear all
set maxvar 32767

/* the data files used in this code include: 

a. randhrs1992_2020v2 can be download from https://hrsdata.isr.umich.edu/data-products/rand-hrs-longitudinal-file-2022

b. A tracker file: trk2018tr_r (https://hrsdata.isr.umich.edu/data-products/cross-wave-tracker-file)

c. Polygenic score files: PGENSCORE3A_R and PGENSCORE3E_R (https://hrs.isr.umich.edu/data-products/genetic-data/products#gdv3)

d. Cross-wave childhood health and family aggregated data: AGGCHLDFH2016A_R
 (https://hrsdata.isr.umich.edu/data-products/cross-wave-childhood-health-and-family-aggregated-data?_gl=1*1e0z3xf*_ga*OTAyMzIwMTI1LjE3NDIzNjkxNjE.*_ga_FF28MW3MW2*czE3NDk1MjYxNjEkbzUkZzEkdDE3NDk1MjY1NTIkajM5JGwwJGgw)
*/ 


global data /Users/yiyue/Library/CloudStorage/OneDrive-CUHK-Shenzhen/data/HRS //xxxx is the path where RAND 1992-2020 longitudinal data file is stored

global working /Users/yiyue/Library/CloudStorage/OneDrive-CUHK-Shenzhen/data/HRS/working

use "$data/randhrs1992_2020v2_STATA/randhrs1992_2020v2",clear

 
merge 1:1 hhid pn using "$data/trk2018v2a/trk2018tr_r"
 drop if _merge <3 
  drop _merge

  g HHID = hhid 
  g PN = pn
  
   
  merge 1:1 HHID PN using "$data/PGENSCORE3/data/PGENSCORE3A_R" // merge with PGS of European ancestry
    drop _merge 
	
   merge 1:1 HHID PN using "$data/PGENSCORE3/data/PGENSCORE3E_R" , update  // merge with PGS of African ancestry 
   drop _merge
   
  
	
	 /* pgs*/ 
	
  drop *PGC10
   
  rename *PGS3_*_* *PGS3_*
  rename AA_PGS3_* *
  
  foreach v of varlist (GENCOG - ANXCC) {
     replace `v' = EA_PGS3_`v' if `v' ==. 
	 }
	 
	 drop EA_PGS3_*
	 
	foreach v of varlist (GENCOG - ANXCC)  {
	  rename `v' PGS_`v'
	  }
		 
	
 
  
/*** the data need to be shaped in long
 format */ 
 keep  r*bmi PGS*    r*diab birthyr ///
   race gender raeduc PC*_* hhid hhidpn   HHID PN  radyear
	  
  
  rename r?bmi rbmi?
  rename r??bmi rbmi??
  
  rename r?diab rdiab?
  rename r??diab rdiab??
  
  
reshape long rbmi rdiab , ///
	i(HHID PN   birthyr  ) ///
	j(wave)
	 
	 
	 g year = 1992+(wave-1)*2
	
	save "$working/HRS92_20",replace
	
	
		 
		/***** cross-wave childhood health and family aggregated data*****/ 
		use "$data/ChildhoodHealthAndFamily/AGGCHLDFH2016A_R",clear
		 foreach v of varlist (FAMFIN- CHOTHCON) {
		  rename `v' ec_`v'
		  }
		  
		  merge 1:m HHID PN using "$working/HRS92_20"
		      drop _merge 
			  
			  	  bys hhidpn wave: drop if _n >1
			  
			 save "$working/HRS92_20",replace 
			 
			 
/****************************************************************************
        variable generating and coding 
****************************************************************************/
		
	/***** construct early life indicator ******/
	//self-rated health 
	recode ec_RTHLTHCH (1=5) (2=4) (3=3) (4=2) (5=1) (8 9=.) 
	  label define healthch 1"poor" 2"fair" 3"above average" 4"very good" 5"excellent" ,modify
	 label values ec_RTHLTHCH healthch
	
	// self-rated family situation
	recode ec_FAMFIN (1=5)  (3=3) (6=2) (5=1) (8 9=.) 
	  label define famfin 1"poor" 2"varied" 3"about average" 5"well off" ,modify
	  label values ec_FAMFIN famfin
	
	//father occupation    ******************** check 
	recode ec_FJOB (1=5) (2=4) (3 6=3) (4=2) (5=1) (8 9 =.) 
	label define fjob 1"manual/operators" 2"service" 3"clerical or armed forces" 4"sales" 5"managerial/professional" ,modify
	   label values  ec_FJOB fjob 
	   
	 // father education 
	 recode ec_FAEDUC (97 98 99=.) 
	 
	 // mother education 
	 recode ec_MOEDUC (97 98 99=.) 
	 
	 // move due to financial difficulty 
	 recode ec_MOVFIN (1=0) (5=1) (8 9=.) 
	 
	 // father unemployed    ************** check 
	  recode ec_FAUNEM (1=1) (5=2) (6 7=0) ( 8 9=.) 
	  
	  // live with grand parents   ********** check 
	  recode ec_LIVEGPAR (1=0) (5=1) (8 9 =.) 
	  
	  // good relationship with father 
	  recode ec_RELWFA (6=.) 
	  
	  // family problems caused by parental drinking 
	  recode ec_DRKDRUG (1=0) (5=1) 
	  
	  // physically abused 
	  recode ec_PHYABUSE (1=0) (5=1) 
	  
	g ec2 = ec_RTHLTHCH+ec_FAMFIN+ec_FJOB+ec_FAEDUC+ec_MOEDUC + ///
	  ec_MOVFIN+ec_FAUNEM+ec_LIVEGPAR+ec_RELWFA +ec_DRKDRUG+ec_PHYABUSE
	  
	  label var ec2 "augmented early life indicator" 
	  
	   save "$working/HRS92_20",replace 
	  
	 
  egen ec2_std = std(ec2)
  g age = year - birthyr
  
	sort hhidpn wave
	 by  hhidpn:  g age_base = age if _n ==1
	  bys hhidpn: egen age_base_2 = mean(age_base)
	  drop age_base
	  rename age_base_2 age_base
	  

	
	 sort hhidpn wave 
	  bys hhidpn: g wave_n = _n
	  
	  
	 mark dead if year >=radyear  | year +2>= radyear 
	 
	 g age_45 = age-45
     
   recode raeduc (1 2=0) (3=1),gen(ba)
   
   recode rbmi (min/29.9=0) (30/max=1), gen(robe)
	  
	  recode rdiab (3=1) (4=0) (5 6=.)
	  
	  g rdiab_raw =rdiab 
	  bys hhidpn: replace rdiab =  1 if (rdiab==0 & rdiab[_n-1] ==1 )  |  (rdiab==0 & rdiab[_n-1]==. & rdiab[_n-2] ==1 )
	   
	   
drop if ec2==.
drop if PGS_BMI ==. 
drop if PGS_T2D ==. 
keep if race==1
	  
	 save  "$working/sample_t2dobe",replace
	 
	 
	 
