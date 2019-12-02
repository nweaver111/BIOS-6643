proc import datafile="'/folders/myfolders/Project/Theoph.csv"
     out=longf/* name for data set for SAS to reference */
     dbms=csv /* identify file as csv */
     replace; /* overwrite HW1_chol if already present */
     getnames=yes; /* take first row as column names from data */
RUN;

ods graphics on;

/* See spaghetti plot of the data: */
PROC SGplot data = longf;
	scatter x=Time y=conc;
	series x=Time y=conc / group=Subject;
RUN;

/* Create 4 columns that will be used for the piecewise cubic regression:
	time_bar: time value providing the maximum concentration for each subject
	s1-s3: cubic components for regression after 3 hours (subjectively determined by looking at spaghetti plots)
*/
DATA longg;
	set longf;
	time_bar = max(IdenticalTime-Time_star, 0);
	if IdenticalTime<25;
	s1 = (max(0,IdenticalTime-3))**1;
	s2 = (max(0,IdenticalTime-3))**2;
	s3 = (max(0,IdenticalTime-3))**3;
RUN;
PROC PRINT data = longg;

/* Run a piecewise cubic regression using subjectively determined knot: */
PROC MIXED data = longg plots=all;
	class Subject;
	model conc = IdenticalTime IdenticalTime*IdenticalTime IdenticalTime*IdenticalTime*IdenticalTime s1 s2 s3 Dose Wt 
	Dose*IdenticalTime Dose*IdenticalTime*IdenticalTime Dose*IdenticalTime*IdenticalTime*IdenticalTime/ solution outp=outer;
	random IdenticalTime/ solution subject = Subject;
RUN;

/* Run the model with maximum for each observation as the knot: */
PROC MIXED data = longg plots=all;
	class Subject;
	model conc = IdenticalTime IdenticalTime*IdenticalTime IdenticalTime*IdenticalTime*IdenticalTime time_bar time_bar*time_bar time_bar*time_bar*time_bar Dose Wt/ solution outp=outer;
	random IdenticalTime/ solution subject = Subject;
RUN;

/* Run the model with linearized gamma as the function */
DATA longg2;
	set longg;
	if IdenticalTime = 0 then delete; /*model wouldn't run with time = 0*/
RUN;

proc nlmixed data=longg2 gconv=1e-10; 
	parms res=2 a1 = 2 a2=.5 a3=9 g1=2; /*Determined by plotting the function on Desmos and messing with parameters until the curve matched my spaghetti plots */
	n=(a3+b1)*(IdenticalTime**(a1))*(exp(-IdenticalTime*(a2))) + a4*Dose + a5*Wt + b2; 
	model conc~normal(n,res); 
	RANDOM b1 b2 ~ NORMAL([0,0],[g1,0,g2]) subject= Subject OUT=re; /*Assume RE's have a simple structure*/
run;

/*run the model with pharmacokinetic one-compartment model as described by Davidian and Giltman: */

PROC NLMIXED data = longg gconv=1e-10;
	PARMS a1 =.6 a2=.1 a3=15; /*determined by plotting the function on Desmos and messing with parameters until the curve matched my spaghetti plots */
	absorp = a1 + b1;
	excret = a2 + b2;
	rati = a3 + b3;
	n=(absorp)*Dose*rati*(absorp-excret)**(-1)*(exp(-excret*IdenticalTime)-exp(-absorp*IdenticalTime));
	MODEL conc~NORMAL(n,res);
	RANDOM b1 b2 b3 ~ NORMAL([0,0,0], [g1, 0, g2, 0, 0, g3]) subject = Subject OUT=re2;
RUN;