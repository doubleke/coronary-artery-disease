libname library '\\tsclient\MacSSD\Users\doubleke\Desktop\PM511';
data cleveland;set library.cad_cleveland_322;run;
data longbeach;set library.cad_longbeach_322;run;
data budapest;set library.cad_budapest_322;run;

proc format;
value malef 1='male' 0='female';
value cpf 1='Typical angina' 2='Atypical angina' 3='Non-anginal pain' 4='Asymptomatic';
value fbsf 0='<120 mg/dl' 1='>=120 mg/dl';
value restecgf 0='Normal' 1='ST-T wave abnormality' 2='Left ventricular hypertrophy';
value exangf 1='yes' 0='no';
value numf 0='No arteries' 1='1 artery' 2='2 arteries' 3='3 arteries' 4='4+ arteries';
run;

data all;
set cleveland longbeach budapest;
label cp='Chest pain'  sysbp='Resting systolic blood pressure, mmHg' fbs='Fasting blood sugar'  restecg='Resting electrocardiographic results' maxhr='Maximum heart rate (bpm)' exang='Exercise induced angina' num='Coronary artery disease (number of arteries with > 50% diameter narrowing)';
format male malef. cp cpf. fbs fbsf. restecg restecgf. exang exangf. num numf.; 
by id;
tempid=1;
run;
/*proc print data=all;run;*/
proc sort data=all; by site;

/*1A*/
data all; set all;
if cp ne . then do; if cp=4 then cpbin=1; else cpbin=0; end;
if restecg ne . then do; if restecg=0 then restecgbin=1; else restecgbin=0; end;
if num ne . then do; if num=0 then numbin=0; else numbin=1; end;
run;

%macro means1(Y);
	proc means data=all n prt std mean missing; 
	class site; 
	var &Y; run;
%mend means1;

%macro anova1(Y);
	proc glm data=all;
	class site;
	model &Y = site/solution; run;

%mend anova1;

title ' age - Mean, n and stdev';
%means1(age);
title ' age - comparison among sites';
%anova1(age);

title ' sysbp - Mean, n and stdev';
%means1(sysbp);
title ' sysbp - comparison among sites';
%anova1(sysbp);

title ' maxhr - Mean, n and stdev';
%means1(maxhr);
title ' maxhr - comparison among sites';
%anova1(maxhr);

proc freq data=all;tables male*site/chisq measures missing;format male malef.; 
title'Male - Freq, n, and comparison among sites';run;


proc freq data=all;tables fbs*site/chisq measures missing;format fbs fbsf.;
title'Fbs - Freq, n, and comparison among sites';run;


proc freq data=all;tables exang*site/chisq measures missing;format exang exangf.;
title'exang - Freq, n, and comparison among sites';run;


proc freq data=all;tables cp*site/missing;format cp cpf.;run;
proc freq data=all;tables cpbin*site/chisq measures missing;format cp cpf.;
title'cp(asymptotic) - Freq, n, and comparison among sites';run;


proc freq data=all;tables restecg*site/missing;format restecg restecgf.;run;
proc freq data=all;tables restecgbin*site/chisq measures missing;format restecg restecgf.;
title'restecgbin(normal) - Freq, n, and comparison among sites';run;

proc freq data=all;tables numbin*site/missing;format num numf.;run;
proc freq data=all;tables numbin*site/chisq measures missing;format num numf.;
title'num(>=1 arteries) - Freq, n, and comparison among sites';run;


/*2A*/
/*Check correlation*/
proc corr data=cleveland plots=matrix;
var maxhr;
with age;
title 'correlation of maxhr and age';
run;
data clevelandnew; set cleveland; age2=age**2; run;
/*Line and outliers and LINE assumptions*/
proc reg data=clevelandnew;
model maxhr=age/r influence;
title 'Regression of Maxhr on age';
output out=regout r=resid p=pred student=student rstudent=jackknife;
run;
proc gplot data=regout; symbol1 v=star; plot (resid student jackknife)*pred/vref=0;run;
title'Assumptions check: maxhr=age';
proc univariate data=regout plot normal; title'Normality';var resid;run;
proc sgplot data=regout; scatter y=resid x=pred; refline 0; loess x=pred y=resid/ clm smooth=.4;title'Homoskedesticity';run;
proc sgplot data=regout; reg y=resid x=maxhr; refline 0; title'Linearity';run;

/*Line and outliers with age and age2, and LINE assumptions*/
proc reg data=clevelandnew;
model maxhr=age age2/r influence;
title 'Regression of Maxhr on age, age2';
output out=regout2 r=resid p=pred student=student rstudent=jackknife;
run;
proc gplot data=regout2; symbol1 v=star; plot (resid student jackknife)*pred/vref=0;run;
title'Assumptions check: maxhr=age';
proc univariate data=regout2 plot normal; var resid;title'Normality';run;
proc sgplot data=regout2; scatter y=resid x=pred; refline 0; loess x=pred y=resid/ clm smooth=.4;title'Homoskedesticity';run;
proc sgplot data=regout2; reg y=resid x=maxhr; refline 0; title'Linearity';run;

/*fit final model*/
proc reg data=clevelandnew outest=trainreg;
maxhr_hat: model maxhr=age age2;title'fit final model';run;

/*2B*/
data longbeachnew; set longbeach; age2=age**2; run;
/*Calculate Y_hat*/
proc score data=longbeachnew score=trainreg out=validreg type=parms nostd predict;
var age age2;run;
/*corr of y and y_hat*/
proc corr data=validreg outp=corrlongbeach; var maxhr maxhr_hat;run;
/*proc print data=corrlongbeach; run;*/
data a; set corrlongbeach; if _TYPE_='CORR' and _NAME_='maxhr'; rename maxhr_hat=PearsonR;run;
/*proc print data=a;run;*/
data b; set a; keep PearsonR R2; R2=PearsonR**2;title'Validation of Long Beach: Correlation between Maxhr and Maxhr_hat';run;
proc print data=b;run;

/*2C*/
/*center age and sysbp*/
data cleveland1; set cleveland; junkid=1; run;
proc means data=cleveland1 noprint; by junkid; var age; output out=outage mean=meanage;
proc means data=cleveland1 noprint; by junkid; var age; output out=outsysbp mean=meansysbp;
data cleveland2; merge outage outsysbp cleveland1; by junkid; 
	agecenter=age-meanage; agecenter2=agecenter**2;
	sysbpcenter=sysbp-meansysbp;sysbpcenter2=sysbpcenter**2;run;

/*Model selection: step 1 - Correlation table*/
proc corr data=cleveland2 best=10; Var agecenter agecenter2 sysbpcenter sysbpcenter2 cp male fbs restecg exang num; 
title'Model selection: step 1 - Correlation table: Correlation and Collinearity';run;

/*Model selection: step 2 - potential interaction terms*/

%macro interaction1(Y);
	proc glm data=cleveland2; model maxhr=agecenter agecenter2 &Y agecenter*&Y /solution; run;
%mend interaction1;

%interaction1(cp);*0.2501;%interaction1(male);*0.5355;%interaction1(sysbpcenter);*0.7680;
%interaction1(fbs);*0.9736;%interaction1(restecg);*0.4356;%interaction1(exang);*0.0199;%interaction1(num);*0.0029;
run;

%macro interaction2(Y);
	proc glm data=cleveland2; model maxhr=agecenter agecenter2 sysbpcenter &Y sysbpcenter*&Y /solution; run;
%mend interaction2;

%interaction2(cp);*0.8525;%interaction2(male);*0.5112;%interaction2(fbs);*0.3507;
%interaction2(restecg);*0.7989;%interaction2(exang);*0.9045;%interaction2(num);*0.8436;
run;


/*Model selection: Step 3 - LINE assumptions and outliers of max model*/
proc glm data=cleveland2;class cp restecg num; 
model maxhr=agecenter agecenter2 agecenter*exang agecenter*num cp male sysbpcenter sysbpcenter2 fbs restecg exang num/solution;
output out=regout r=resid p=pred student=student rstudent=jackknife;
run;
Title 'Maximun Model';
proc gplot data=regout; symbol1 v=star; plot (resid student jackknife)*pred/vref=0;run;

proc univariate data=regout plot normal; var resid;title'Normality';run;
proc sgplot data=regout; scatter y=resid x=pred; refline 0; loess x=pred y=resid/ clm smooth=.4;title'Homoskedesticity';run;
proc sgplot data=regout; reg y=resid x=maxhr; refline 0; title'Linearity';run;

/*Model selection: Step 4 - Stepwise method*/
proc glmselect data=cleveland2; class cp restecg num; 
model maxhr=agecenter agecenter2 agecenter*exang agecenter*num cp male sysbpcenter sysbpcenter2 fbs restecg exang num/selection=stepwise details=all hierarchy=single sls=0.15 sle=0.15 showpvalues select=sl stop=sl;
run;

/*Model selection: Step 5 - report and check final model */
proc glm data=cleveland2;
model maxhr=agecenter agecenter2 agecenter*exang sysbpcenter exang num/ solution;
output out=j1 r=resid p=pred student=student rstudent=jackknife;
title'Final model';
run;*r2=0.6533;
proc gplot data=j1; symbol1 v=star; plot (resid student jackknife)*pred/vref=0;run;
proc univariate data=j1 plot normal; var resid;title'Normality';run;
proc sgplot data=j1; scatter y=resid x=pred; refline 0; loess x=pred y=resid/ clm smooth=.4;title'Homoskedesticity';run;
proc sgplot data=j1; reg y=resid x=maxhr; refline 0; title'Linearity';run; 

/*2D*/
/*Validation*/
data cleveland2; set cleveland2; agecenterexang=agecenter*exang; run;

/*fit final model*/
proc reg data=cleveland2 outest=trainreg;
maxhr_hat: model maxhr=agecenter agecenter2 agecenterexang sysbpcenter exang num;
title 'Fit final model';run;

%macro validation(dataset);
/*Validation with newdataset*/
/*center age and sysbp*/
data temp1; set &dataset; junkid=1; run;
proc means data=temp1 noprint; by junkid; var age; output out=outage mean=meanage;
proc means data=temp1 noprint; by junkid; var sysbp; output out=outsysbp mean=meansysbp;
data temp2; merge outage outsysbp temp1; by junkid; 
	agecenter=age-meanage; agecenter2=agecenter**2;
	sysbpcenter=sysbp-meansysbp;
	agecenterexang=agecenter*exang;
run;

/*Calculate Y_hat*/
proc score data=temp2 score=trainreg out=validreg type=parms nostd predict;
var agecenter agecenter2 agecenterexang sysbpcenter exang num;run;

/*corr of y and y_hat*/
proc corr data=validreg outp=corrvalid; var maxhr maxhr_hat;run;
/*proc print data=corrvalid; run;*/
data a; set corrvalid; if _TYPE_='CORR' and _NAME_='maxhr'; rename maxhr_hat=PearsonR;run;
/*proc print data=a;run;*/
data b; set a; keep PearsonR R2; R2=PearsonR**2;run;
proc print data=b; run;

%mend validation;

title 'Validation: R2 with dataset = Long Beach';
%validation(longbeach);
title 'Validation: R2 with dataset = Budapest';
%validation(budapest);



/*3*/
proc means data=all noprint; by tempid; var age; output out=outage mean=meanage;
data all2; merge outage all; by tempid; agecenter=age-meanage;run;
proc sort data=all2; by site;run;

data num0; set all2; if num=0; run;
data num1; set all2; if num>=1; run;


/*Confounder1? Agecenter*/
%macro confound1(dataset);
	proc corr data=&dataset; var sysbp agecenter; title 'Confounder1: corr between Sysbp and potential confounder: age'; run;/*Y and Confounder*/

	proc glm data=&dataset;
	class male;
	model agecenter= male/solution;title 'Confounder1: corr between gender and potential confounder: age';  run;/*X and confounder*/
%mend confound1;

/*Confounder2? Site*/
%macro confound2(dataset);
	proc glm data=&dataset;
	class site;
	model sysbp=site;title 'Confounder2: corr between Sysbp and potential confounder: site'; run;/*Y and Confounder*/

	proc freq data=&dataset;
	tables site*male/chisq measures missing cmh;
	title 'Confounder2: corr between gender and potential confounder: site';/*X and confounder*/
	run;
%mend confound2;


/*3A*/
%confound1(num0);%confound2(num0);
proc glm data=num0; 
class male; model sysbp=male/solution;
title 'ANOVA: Sysbp difference between genders'; 
title2'patients who do not have coronary artery disease (num=0)';
run;

/*3B*/
%confound1(num1);%confound2(num1);
proc glm data=num0; 
class male; model sysbp=male/solution;
title 'ANOVA: Sysbp difference between genders'; 
title2'patients who have coronary artery disease (num>=1)';
run;



/*4A*/

/*population of interest*/
data eligible; set budapest; if sysbp ne . and sysbp>140; if male=1; run;
/*proc print data=eligible;title 'population of interest'; run;*/
/*Step 1 - get sample of 20*/
proc surveyselect data=eligible out=sample20 method=SRS seed=8186 sampsize=20;
id _all_;title 'Step 1 - sample of 20 men, with high sysbp from Budapest';
run;
proc print data=sample20;title 'sample of 20 men, with high sysbp from Budapest'; run;

proc sort data=sample20; by id;run;

/*Step 2 - Grouping*/
proc surveyselect data=sample20 out=group1 method=SRS seed=2016 sampsize=10;
id _all_;title 'Step 2 - Assign treatment group';
run;
data treatment; set group1; treatment=1; run;
/*proc print data=treatment;title 'treatment group'; run;*/

proc sort data=treatment; by id;run;

data placebo; merge sample20(in=inall) treatment(in=intreat);
by id; if intreat=0 and inall=1; treatment=0; run;
/*proc print data=placebo; title 'placebo group';run;*/

/*Step 3 - summary*/
data recruited(keep=id age male sysbp treatment); merge placebo treatment; by id; run;
proc print data=recruited; title 'Step 3 - result: recruited participants and treatment assignments'; run;


/*4B*/
proc means data=eligible n std mean; var sysbp; title'Get stdev. of the population of interest';run;

data spower; 
delta=10; sigma=12.268; /*delta is the diffrence of mean; sigma is std of eligible men from Budapest*/
alpha=0.05; sides=2; N=20; 
Qe=0.5; Qc=1-Qe; /*Qe is the proportion of treatment group*/ 
Zalpha=probit(1-alpha/sides);
Zbeta=sqrt(N*delta**2/(sigma**2*(1/Qe+1/Qc)))-Zalpha;
power=probnorm(Zbeta);
run;
proc print data=spower; title 'Sample power';run;
