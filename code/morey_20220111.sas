
libname eff "U:\Consulting\KEL\Leach\Morey\data"; 
%include "U:\Consulting\KEL\Fox\Weisse\Cote\code\ROC_Optimal_Cutoff_031816.sas";

data eff /* (drop=itapx) */; set eff.hemo_effusion_20220114 /* (rename=iTAPSE=itapx) */; 
c_nc = cardiac_versus_noncardiacdx;
if cardiac_versus_noncardiacdx ne . then do;
c_nc2 = (cardiac_versus_noncardiacdx in (1 2));
end;

/*iTAPSE = input(itapx,5.0);*/
run;


data eff23; set eff;
where c_nc ne 1;
run;

/*proc freq;*/
/*tables cardiac_versus_noncardiacdx*(c_nc c_nc2);*/
/*run;*/

%let contvars = HCT WBC PLT BUN Creat Tot_pro Alb Glob ALT ALKP GGT Tbili iTAPSE FAC RAD NT_proBNP cTnI;
%let catvars = Heart_murmur Jug_dist Jug_puls Abd_dist Fluid_wave 
Hep_venous_dist_RV1 Abd_eff_RV1 CaudalVC_dist_RV1 Peri_eff_RV1 GBW_oedema_RV1 Pl_eff_RV1 
Hep_venous_dist_RV2 Abd_eff_RV2 CaudalVC_dist_RV2 Peri_eff_RV2 GBW_oedema_RV2 Pl_eff_RV2;
/*Hep_venous_dist Abd_eff CaudalVC_dist Peri_eff GBW_oedema Pl_eff;*/
/*NT_proBNP and cTnI*/

%let ds = eff;

/* Descriptive statistics by group (Table 1) Continuous with p-values */
%macro desc_cont(iv,model,research_q);
data ds; set &ds.; run;

%let vars = &contvars.;
%let word_cnt=%sysfunc(countw(&vars));
%do i = 1 %to &word_cnt.;
/*%do i = 1 %to 1;*/


proc means data = ds n nmiss min max mean p25 median p75 std var ;
class &iv; 
var %scan(&vars.,&i);
ods output summary = xxx;
run;
data m_%scan(&vars.,&i);
length dv iv $50.;
set xxx;
rename %scan(&vars.,&i)_n = n;
rename %scan(&vars.,&i)_NMiss = n_Miss;
rename %scan(&vars.,&i)_mean = mean;
rename %scan(&vars.,&i)_stddev = stddev;
rename %scan(&vars.,&i)_var = variance;
rename %scan(&vars.,&i)_min = min;
rename %scan(&vars.,&i)_max = max;
rename %scan(&vars.,&i)_p25 = p25;
rename %scan(&vars.,&i)_median = median;
rename %scan(&vars.,&i)_p75 = p75;
dv = "%scan(&vars.,&i)"; iv = "&iv."; research_q = "&research_q."; model = "&model."; 
run;
data rawmean;
set rawmean m_%scan(&vars.,&i) (rename=&iv.=iv_level);
run;
proc datasets lib=work nodetails nolist;
delete mean_%scan(&vars.,&i);
run;
proc mixed data=ds;
class &iv. ;
model %scan(&vars.,&i) = &iv. ;
lsmeans &iv. / pdiff =all; 
ods output tests3 = t lsmeans = lsm diffs=d;
data t (rename=probf = p_value); length var $50.;
set  t (keep=ProbF); var = "%scan(&vars.,&i)"; research_q = "&research_q."; model = "&model."; run;
data ls (rename=(estimate = lsmean &iv.=iv_level)); length iv var $50.; iv = "&iv.";
set  lsm (keep=estimate stderr &iv.); var = "%scan(&vars.,&i)"; research_q = "&research_q."; model = "&model."; run;
data d (drop=df tvalue) ; length dv iv $50.;
set  d (rename=(&iv=iv_grp1 _&iv. = iv_grp2 estimate=grp_mn_diff effect=iv)); dv = "%scan(&vars.,&i)"; research_q = "&research_q."; model = "&model."; run;
data t3 ; set t3 t ; data lsmean ; set lsmean ls ; data diffs; set diffs d ;   run;
%end;
%mend;

data rawmean; set _null_; data t3; set _null_; data lsmean; set _null_; data diffs; set _null_; run;
%desc_cont(c_nc,ANOVA,1);
%desc_cont(c_nc2,ANOVA,2);



%macro penorm (ds,iv);
%let vars = &contvars ;
%let count=%sysfunc(countw(&vars.));
%do i = 1 %to &count.;
proc reg data=&ds.;
  model %scan(&vars,&i.)= &iv. ;
  output out=%scan(&vars,&i.)res rstudent=%scan(&vars,&i.)r /*h=lev cookd=cd dffits=dffit*/ ;
  ods output ParameterEstimates=%scan(&vars,&i.)pe;
run;
quit;
data %scan(&vars,&i.)pe;
length var $50;
set %scan(&vars,&i.)pe;
var = "%scan(&vars,&i.)";
where variable not in ("Intercept");
run;
data pe;
set pe %scan(&vars,&i.)pe;
run;
* Examine residuals for normality ;
proc univariate data=%scan(&vars,&i.)res plots plotsize=30 normal;
  var %scan(&vars,&i.)r;
  ods output TestsForNormality=%scan(&vars,&i.)norm;
run;
data %scan(&vars,&i.)norm;
length var $50;
set %scan(&vars,&i.)norm;
var = "%scan(&vars,&i.)";
run;
data norm;
set norm %scan(&vars,&i.)norm;
run;

* Non Parametric;

 proc npar1way data = &ds. median wilcoxon;
 class &iv.;
 var %scan(&vars,&i.);
 ods output mediantest = %scan(&vars,&i.)mt wilcoxontest=%scan(&vars,&i.)wt KruskalWallisTest=%scan(&vars,&i.)kw;
 run;
 data %scan(&vars,&i.)mt;
 length Name1 Label1 cValue1 var $50;
set %scan(&vars,&i.)mt;
var = "%scan(&vars,&i.)";
where name1 in ("PR_MED" "P2_MED" "PL_MED");
run;
 data %scan(&vars,&i.)wt;
 length Name1 Label1 cValue1 var $50;
set %scan(&vars,&i.)wt;
var = "%scan(&vars,&i.)";
where name1 in ("PL_WIL" "P2_WIL" "PR_WIL");
run;
 data %scan(&vars,&i.)kw;
 length Name1 Label1 cValue1 var $50;
set %scan(&vars,&i.)kw;
var = "%scan(&vars,&i.)";
where name1 in ("P_KW");
run;
data npar;
set npar %scan(&vars,&i.)mt %scan(&vars,&i.)wt %scan(&vars,&i.)kw;
run;
	%end;
	data npar_wil; set npar; where name1 = "P2_WIL"; run;
	data norm_ks; set norm; where testlab = "D"; run;
%mend;
data norm; set _null_; run; data pe; set _null_; data pe_norm; set _null_; data npar; set _null_; run;
%let ds = eff;
/*%penorm(&ds.,c_nc);*/
%penorm(&ds.,c_nc2);


%macro printnorm(ds);
proc print data= &ds.; run;
%mend;
title "Normal Test Stats as guide";
%printnorm(norm_ks);
title "Non-parametric output as mean alternative";
%printnorm(npar_wil);


	*** Chi Square Cross-Tab Frequency Anlaysis Set-up;
%macro chi(iv1,model,research_q);

data ds; set &ds.; run;

%let vars = &catvars.;
%let word_cnt=%sysfunc(countw(&vars));
%do i = 1 %to &word_cnt.;
/*%do i = 1 %to 1;*/


proc freq data =  ds;
tables &iv1.*%scan(&vars.,&i) / chisq;
ods output CrossTabFreqs = freqs ChiSq= chi FishersExact = fish ;run;
run;
data freqs_in; set freqs (drop =  _table_ rowpercent colpercent); 
where _type_ in ("11" "00"); drop _type_; run;

data ctf_&research_q. (drop =  &iv1. %scan(&vars.,&i)); retain model /* analysis */ research_q; length model research_q $50.; set freqs_in; 
/* analysis = "&analysis."; */ research_q = "&research_q."; model = "&model."; 
iv1 = put(&iv1.,5.0); iv2 = put(%scan(&vars.,&i),5.0);
run;
data ctf_all; set ctf_all ctf_&research_q.; run;

	%if %sysfunc(exist(fish)) %then %do;
data f_&research_q. (rename=nvalue1=P_FISH); set fish (keep =  table name1 nvalue1); 
where name1 in ("P_TABLE"); drop name1; run;
data f_p; set f_p f_&research_q.; run;
	%end;
	%if %sysfunc(exist(chi)) %then %do;
data chi_&research_q. (rename=prob=P_CHI); length model research_q $50.;
set chi (keep =  table statistic prob); 
where statistic in ("Chi-Square");  drop statistic; 
model = "&model."; research_q = "&research_q."; run;
data chi_p; set chi_p chi_&research_q.; run;
	%end;

proc datasets library=work nolist nodetails; delete chi fish freq p_&research_q. chi_&research_q.; 
run;

%end;
	%mend;

data ctf_all; set _null_; data f_p; set _null_; data chi_p; set _null_; run;

%chi(c_nc,Chi_Square,3);
%chi(c_nc2,Chi_Square,4);

/* Part C) Additional Request for specific Sensitivity / Specificity at 0.9*/

%macro logit(ds,iv,dv,prob);
	/* ods graphics on;*/



	proc logistic data=&ds. /* plots=roc */;
		model &iv. (event = "1") = &dv. / outroc=roc1;
		output out=out p=phat;
	run;

	/*ods graphics off;*/
	title  "ROC plot for &dv. = &iv.";
	title2 " ";

	%rocplot( inroc = roc1, inpred = out, p = phat,
		id = &dv. _sens_ _spec_ _OPTEFF_ _opty_ _cutpt_,
		optcrit = youden , pevent = &prob.,
		optsymbolstyle = size=0, optbyx = panelall, x = &dv.)

		data rocme (drop =  &dv.);
	retain _id &dv. _CORRECT_ _sens_ _spec_ _FALPOS_ _FALNEG_ _opty_;
	set _rocplot (keep=_id &dv. _CORRECT_ _sens_ _spec_    _sensit_ _FALPOS_ _FALNEG_ _opty_ __spec_ _POS_ _NEG_);
dv = "&dv.";
iv = "&iv.";
/*	if (_sens_ in (0 1)) or (_spec_ in (0 1)) or (_opty_ = "Y") then*/
		if (_opty_ = "Y") then
		output;
	run;
/**/
/*	data rocme9_10;*/
/*	retain _id &iv. _CORRECT_ _sens_ _spec_ _FALPOS_ _FALNEG_ _opty_;*/
/*	set _rocplot (keep=_id &dv. _CORRECT_ _sens_ _spec_    _sensit_ _FALPOS_ _FALNEG_ _opty_ __spec_ _POS_ _NEG_);*/
/**/
/*	if (0.89 le _sens_ le 0.91) or (0.09 le _spec_ le 0.11) or (_opty_ = "Y") then*/
/*		output;*/
/*		run;*/
/**/
/*		data rocme90;*/
/*	retain _id &dv. _CORRECT_ _sens_ _spec_ _FALPOS_ _FALNEG_ _opty_;*/
/*	set _rocplot (keep=_id &dv. _CORRECT_ _sens_ _spec_    _sensit_ _FALPOS_ _FALNEG_ _opty_ __spec_ _POS_ _NEG_);*/
/**/
/*	if (0.89 le _sens_ le 0.91) or (0.89 le _spec_ le 0.91) or (_opty_ = "Y") then*/
/*		output;*/
/*		run;*/
data rocme_all; set rocme_all rocme; run;

%mend;


/*%macro logme;*/
/*%let vars = HCT WBC PLT BUN Creat Tot_pro Alb Glob ALT ALKP GGT Tbili iTAPSE FAC RAD NT_proBNP cTnI;*/
/*%let word_cnt=%sysfunc(countw(&vars));*/
/*%do i = 1 %to &word_cnt.;*/
/*/*%do i = 1 %to 1;*/*/
/*%logit(eff,c_nc2,%scan(&vars.,&i),1);*/
/*%end;*/
/*%mend;*/
/*%logme;*/
;




data rocme_all; set _null_; run;

%logit(eff,c_nc2,HCT,1);
%logit(eff,c_nc2,WBC,1);
%logit(eff,c_nc2,PLT,1);
%logit(eff,c_nc2,BUN,1);
%logit(eff,c_nc2,Tot_pro,1);
%logit(eff,c_nc2,Alb,1);
%logit(eff,c_nc2,Glob,1);
%logit(eff,c_nc2,ALT,1);
%logit(eff,c_nc2,ALKP,1);
%logit(eff,c_nc2,GGT,1);
%logit(eff,c_nc2,Tbili,1);
%logit(eff,c_nc2,iTAPSE,1);
%logit(eff,c_nc2,FAC,1);
%logit(eff,c_nc2,RAD,1);
%logit(eff,c_nc2,NT_proBNP,1);
%logit(eff,c_nc2,cTnI,1);

/* Reduced Population omitting c_nc = 1 with new data set eff23*/

data rocme_all; set _null_; run;

%logit(eff23,c_nc2,HCT,1);
%logit(eff23,c_nc2,WBC,1);
%logit(eff23,c_nc2,PLT,1);
%logit(eff23,c_nc2,BUN,1);
%logit(eff23,c_nc2,Tot_pro,1);
%logit(eff23,c_nc2,Alb,1);
%logit(eff23,c_nc2,Glob,1);
%logit(eff23,c_nc2,ALT,1);
%logit(eff23,c_nc2,ALKP,1);
%logit(eff23,c_nc2,GGT,1);
%logit(eff23,c_nc2,Tbili,1);
%logit(eff23,c_nc2,iTAPSE,1);
%logit(eff23,c_nc2,FAC,1);
%logit(eff23,c_nc2,RAD,1);
%logit(eff23,c_nc2,NT_proBNP,1);
%logit(eff23,c_nc2,cTnI,1);




**** K Measure **********;
/*Hep_venous_dist_RV1 Abd_eff_RV1 CaudalVC_dist_RV1 Peri_eff_RV1 GBW_oedema_RV1 Pl_eff_RV1 */
/*Hep_venous_dist_RV2 Abd_eff_RV2 CaudalVC_dist_RV2 Peri_eff_RV2 GBW_oedema_RV2 Pl_eff_RV2*/
/*noprint;*/


%macro kappame(iv1,iv2);
/*ods trace on;*/
proc freq data=eff;
tables &iv1.*&iv2. / agree; test wtkap;
ods output CrossTabFreqs=ctf SymmetryTest=st kappa=k wtkappa=wtk;
run;
data ctf (keep= table iv1 iv2 frequency percent missing); length iv1 iv2 $50.; set ctf; iv1 = "&iv1."; iv2 = "&iv2.";  run;
data st (keep=table sym_pvalue); set st; where name1="P_TSYMM"; sym_pvalue=nvalue1; run;
data k (keep=table kappa_value); set k; kappa_value=nvalue1; where label1="Kappa"; run;
data wtk (keep=wt_kappa_value); set wtk; wt_kappa_value=nvalue1; where name1="_WTKAP_"; run;
data st_wtk_k; merge k wtk st; run;
data ctf_all; set ctf_all ctf; run;
data wtk_k_all; set wtk_k_all st_wtk_k; run;

/*ods trace off;*/
%mend;


data ctf_all; set _null_; data wtk_k_all; set _null_; run;
/*https://www.stattutorials.com/SAS/TUTORIAL-KAPPA.htm*/
%kappame(Hep_venous_dist_RV1,Hep_venous_dist_RV2);
/*%kappame(Abd_eff_RV1,Abd_eff_RV2);*/
%kappame(CaudalVC_dist_RV1,CaudalVC_dist_RV2);
%kappame(Peri_eff_RV1,Peri_eff_RV2);
%kappame(GBW_oedema_RV1,GBW_oedema_RV2);
%kappame(Pl_eff_RV1,Pl_eff_RV2);
    

/*%kappame(Abd_eff_RV1,Abd_eff_RV2); * All 1/1;*/
proc freq data=eff;
tables Abd_eff_RV1*Abd_eff_RV2 / agree chisq; test wtkap;
/*ods output CrossTabFreqs=ctf SymmetryTest=st kappa=k wtkappa=wtk;*/
run;











