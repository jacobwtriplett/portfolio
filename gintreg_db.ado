/*This ado file executes non-linear interval regressions where the error term is 
distributed in the GB2 or SGT family tree

Author--Jacob Orchard
v 1.4

******************************

Update for v 1.5 (5/19/2017) by Will Cockriel
1. Added repeat option 

******************************

Update for v 1.6 (1/12/2018) by Bryan Chia
1. Added no constant option
2. Changed the initial option to cases like the no constant option
3. Took out the limits on the number of initial parameters as it complicated stuff with heteroskedasticity 
4. Initially, I was thinking about having users put in only p,q, lambda if its heteroskedastic and 
   mu, sigma, p, q, lambda if it is homoskedastic. 
   I think that is probably too confusing so I standardized it. Decided to just have them put in p, q and lambda, 
   whichever is relevant. I edited the help file too. 
   Decided that if it is heteroskedastic though it is hard to know what sigma is 
   so that was my thoughts for why it would be good for them to just guess p, q and lambda.
   Mu and sigma is based on normal/ln values depending on what family the distribution is from. 
5. Edited some of the mean values that were incorrect 
6. We use an lnsigma now instead of sigma. All intllf files have been edited to exp(sigma) rather than sigma
   This is potentially confusing. Maybe I should change "sigma" to "lnsigma" for all the functions...
7. Edited nortolerance to nonrtolerance 
8. I took out the constant only ml optimization before the full model as there were some convergence issues
   Now, if they put in initial values, it will instead first optimize assuming the simplest distributions
   in each family, namely the lognormal and the normal and use those values as start values instead. 
9. Edited lambda to make it bounded between 0 to 1. Had to make some changes to transform what we call alpha
   back to lambda. 

******************************

Update for v 1.7 (3/27/2018) by Bryan Chia & Jonny Jensen
1. Fixed pdfs to allow better convergence for point estimates with sgt family
2. Point estimates - change sigma to exp(sigma) 
3. Add BIC and AIC 

******************************

Update for v 1.8 (5/22/2018) by Bryan Chia 
1. Changed it such that users can now specify the following distributions as well:
- for the SGT family: GT, ST, GED, SLaplace, T
- for the GB2 family: Br 12, Br3, Gamma 

******************************

Update for v 2.0 (12/13/2021) by Jacob Triplett
1. Fixed misidentified constraints for the Gamma and Weibull distributions
2. Corrected mechanism by which distribution defaults to Normal
3. BIC and AIC values are now returnable with e(BIC) and e(AIC) and are only displayed when option bicaic
   is specified and for gb2-tree distributions and if data is non-grouped
4. Standardized output title formatting
5. Initiated gini option, currently only coded for weibull, br3, br12 but can accomodate more
6. Conventionalized distribution input names to be consistent with what Stata's in-house programs use
   and trimmed superfluous labels (i.e. dropped "ln" when "lnormal" is also accepted and prefererred 
   conventionally, changed "gg" to "ggamma") 
7. Added plot option to easily vizualize pdf-- works for all accepted distributions
8. Users can now specify the following distributions as well: laplace, skewed normal ("snormal")
*/

capture program drop gintreg_db

program gintreg_db, eclass
version 13.0
	if replay() {
		display "Replay not implemented"
	}
	else {
		di as text "NOTICE: this is an unofficial copy of gintreg for debugging purposes only"
		set more off
		syntax varlist(min=2 fv ts)  [aw fw pw iw] [if] [in] ///
		[, DISTribution(string) /// 
		sigma(varlist) ///
		lambda(varlist) ///
		p(varlist) ///
		q(varlist) ///
		b(varlist) ///
		beta(varlist)  ///
		INITial(numlist) ///
		vce(passthru)  ///
		eyx(string) ///
		Het(string) CONSTraints(passthru) DIFficult TECHnique(passthru) ITERate(passthru)  /// 
		nolog TRace GRADient showstep HESSian SHOWTOLerance TOLerance(passthru) NONRTOLerance ///
		LTOLerance(passthru) NRTOLerance(passthru) robust cluster(passthru) repeat(integer 1) NOCONStant ///
		gini bicaic plot(string) ///
		svy SHOWConstonly FREQuency(varlist)] 
		
		*Defines Independent and Dependent Variables
		local depvar1: word 1 of `varlist'
		local depvar2: word 2 of `varlist'
		local tempregs: list varlist - depvar1 
		local regs: list tempregs - depvar2
				
		*Defines variables for other parameters
		if "`sigma;" != ""{
			local sigmavars `sigma'
			}
			
		if "`lambda;" != ""{
			local lambdavars `lambda'
			}
		if "`p;" != ""{
			local pvars `p'
			}
		if "`q;" != ""{
			local qvars `q'
			}
		
		if "`het'" != "" {
             ParseHet `het'
             local hetvar "`r(varlist)'"
             local hetnocns "`r(constant)'"		
			 }
		
		
		local nregs: word count `regs'
		local nsigma: word count `sigmavars'
		local nlambda: word count `lambdavars'
		local np: word count `pvars'
		local nq: word count `qvars'
		
		*Working with heteroskedasticity
		
		if "`het'" != "" {
		local sigmaeq `"(`hetvar')"'
		di as txt "`sigmaeq'"
		}
		
		
		
		*Displays error if distribution name not in library
		if `"`distribution'"' != "" {
			if "`distribution'" != "lnormal" & "`distribution'" != "br12" & "`distribution'" != "br3" & ///
		"`distribution'" != "ggamma" & "`distribution'" != "gb2" & "`distribution'" != "gamma" /// 
		& "`distribution'" != "sgt" & "`distribution'" != "gt" & "`distribution'" != "st" & "`distribution'" != "ged" ///
		& "`distribution'" != "sged" & "`distribution'" != "t" & "`distribution'" != "slaplace" & "`distribution'" != "laplace" & "`distribution'" != "weibull" ///
		& "`distribution'" != "normal" & "`distribution'" != "snormal"  {
			di as err "Distribution specified incorrectly. Use the following: gb2, br12, br3, gamma, ggamma, weibull, sgt, gt, st, sged, ged, laplace, slaplace, t, normal, snormal, lnormal"
			exit 498
			}
		}
		
		
		*Displays error if using the wrong parameter with chosen distribution
		if  (`nlambda' > 0) & (("`distribution'" != "sgt") & ("`distribution'" != "sged")){
				di as err "Lambda is not a parameter of the chosen distribution"  
				exit 498 
			}
			
		if `np' > 0 & ("`distribution'" != "sgt" & "`distribution'" != "gb2" & "`distribution'" != "ggamma" ///
					& "`distribution'" != "sged" & "`distribution'" != "gamma" & "`distribution'" != "br3") {
					di as err "p is not a parameter of the chosen distribution"  
					exit 498
				}
		if  `nq' > 0 &  ("`distribution'" != "sgt" & "`distribution'" != "gb2" & "`distribution'" != "weibull" & "`distribution'" != "br12") {

						di as err "q is not a parameter of the chosen distribution"
						exit 498 
				}
				
				
		*Displays error if depvar1 is greater than depvar2
		qui count if `depvar1' > `depvar2' & `depvar1' != .
		if r(N) >0{
			di as err "Dependent variable 1 is greater than dependent variable 2 for some observation"
			exit 198
		}
		
		*Defines titles used when running the program
	    local gb2title "{title:Interval Regression with GB2 Distribution}"
		local ggtitle "{title:Interval Regression with Generalized Gamma Distribution}"
	    local lntitle "{title:Interval Regression with Log-Normal Distribution}"
		local sgttitle "{title:Interval Regression with SGT Distribution}"
		local sgedtitle "{title:Interval Regression with the SGED Distribution}"
		local slaplacetitle "{title:Interval Regression with the Skewed Laplace Distribution}"
		local normaltitle "{title:Interval Regression with Normal Distribution}"
		local snormaltitle "{title:Interval Regression with Skewed Normal Distribution}"
		local br12title "{title:Interval Regression with Burr 12 Distribution}" 
		local br3title "{title:Interval Regression with Burr 3 Distribution}"
		local gammatitle "{title:Interval Regression with Gamma Distribution}"
		local gttitle "{title:Interval Regression with GT Distribution}"
		local sttitle "{title:Interval Regression with ST Distribution}"
		local gedtitle "{title:Interval Regression with GED Distribution}"
		local ttitle "{title:Interval Regression with T Distribution}"
		local laplacetitle "{title:Interval Regression with Laplace Distribution}"
		local weibulltitle "{title:Interval Regression with Weibull Distribution}" 
		
		*Decides which observations to use in analysis.
		
		marksample touse, nov
		
		foreach i in  `regs' `sigmavars' `pvars' `qvars'{
			qui replace `touse' = 0 if `i' ==.
			}
		qui replace `touse' = 0 if `depvar1' == `depvar2' == .
		
		*Gets rid of uncensored observations with a non-positive dependent
		*variable if user is using a positive distribution.
	
		if "`distribution'" == "lnormal" | "`distribution'" == "br12"| "`distribution'" == "br3"| ///
		"`distribution'" == "ggamma" | "`distribution'" == "gb2" | "`distribution'" == "gamma" | ///
		"`distribution'" == "weibull" {
			quietly{ 
			  count if `depvar1' < 0 & `touse' & `depvar1' == `depvar2'
			  local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				if `n' == 1{
					noi di as txt " {res:`depvar1'} has `n' uncensored value < 0;" _c
					noi di as text " not used in calculations"
				}
				else{
					noi di as txt " {res:`depvar1'} has `n' uncensored values < 0;" _c
					noi di as text " not used in calculations"
					}
				}

			  count if `depvar1' == 0 & `touse' & `depvar1' == `depvar2'
			  local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				noi di as txt " {res:`depvar1'} has `n' uncensored values = 0;" _c
				noi di as text " not used in calculations"
				}
				
			  count if `depvar1' <= 0 & `depvar2' <= 0 & `touse' & `depvar1' != `depvar2'
			  local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				noi di as txt " {res:`depvar1'} has `n' intervals < 0;" _c
				noi di as text " not used in calculations"
				}
				
			count if `depvar1' == . & `depvar2' <= 0 & `touse' & `depvar1' != `depvar2'
			local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				noi di as txt " {res:`depvar1'} has `n' left censored values <= 0;" _c
				noi di as text " not used in calculations"
				}
				
		  replace `touse' = 0 if  `depvar2' <= 0
		  
		  }
		}
		
		*Counts the number of each type of interval
		quietly{
			count
			local total = r(N)
			count if `depvar1' != . & `depvar2' != . & `depvar1' == `depvar2'  /// 
			& `touse' == 1
			local nuncensored = r(N)
			count if `depvar1' != . & `depvar2' != . & `depvar1' != `depvar2'  ///
			& `touse' == 1
			local ninterval = r(N)
			count if `depvar1' != . & `depvar2' == .  & `touse' == 1
			local nright = r(N)
			count if `depvar1' == . & `depvar2' != . & `touse' == 1
			local nleft = r(N)
			count if `depvar1' == . & `depvar2' ==. & `touse' == 1
			local nnoobs = r(N)
			
		}
		*Duplicates observations if group data
		if "`frequency'" != ""{
			tempvar tot per
			qui egen `tot' = sum(`frequency')
			qui gen `per' = `frequency'/`tot'
			global group_per `per'
		}
		
		*Evaluates model if grouped data
if "`frequency'" != ""{	

		if "`distribution'" == "normal" {
		
			if "`noconstant'" != ""{
			
			local evaluator intllf_normal_group
			
			di " "
			di as txt "Fitting Full model with no constant:"
			
			
			if "`initial'" !=""{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial',copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			}
			
			else{
			local evaluator intllf_normal_group
			
			if "`initial'" !=""{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial', copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display			
			}
			}
		}
		
		else if "`distribution'" == "lnormal" {
		
			if "`noconstant'" != ""{
			
			local evaluator intllf_ln_group
			
			di " "
			di as txt "Fitting Full model with no constant:"
			
			if "`initial'" !=""{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial',copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			initial(`initial') `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			}
			
			else{
			local evaluator intllf_ln_group
			
			if "`initial'" !=""{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial', copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display			
			}
			}
		}
		
		else if "`distribution'" == "t" {
		
			constraint define 1 [p]_cons=2
			constraint define 2 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"

			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
		
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "st" {
		
			constraint define 1 [p]_cons=2
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"

			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
		
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "gt" {
		
		constraint define 1 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "sgt" {
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sgttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sgttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "laplace" {
		
		constraint define 1 [p]_cons=1
		constraint define 2 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`laplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`laplacetitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`laplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`laplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "slaplace" {
		
		constraint define 1 [p]_cons=1
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`slaplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`slaplacetitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`slaplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`slaplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "ged" {
		
		constraint define 1 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gedtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "sged" {
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sgedtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "snormal" {
			
			constraint define 1 [p]_cons=2
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`snormaltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`snormaltitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`snormaltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`snormaltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "br12" {
		
		constraint define 1 [p]_cons=1

			if "`noconstant'" != ""{
			
			local evaluator intllf_gb2exp_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(`initial' coeff2,copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br12title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			qui ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br12title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_gb2exp_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br12title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			
			di " "
			di as txt "Fitting Full model:"
						
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br12title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
				}
			}
			
			}
		
		else if "`distribution'" == "br3" {
		
		constraint define 1 [q]_cons=1

			if "`noconstant'" != ""{
			
			local evaluator intllf_gb2exp_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(`initial' coeff2,copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br3title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			qui ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br3title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_gb2exp_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br3title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			
			di " "
			di as txt "Fitting Full model:"
						
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br3title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
				}
			}
			
			}
			
		else if "`distribution'" == "gb2" {

			if "`noconstant'" != ""{
			
			local evaluator intllf_gb2exp_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(`initial' coeff2,copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			qui ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gb2title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_gb2exp_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			
			di " "
			di as txt "Fitting Full model:"
						
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gb2title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
				}
			}
			
			}
		
		else if "`distribution'" == "gamma" {

		constraint define 1 [lnsigma]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ggsigma_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gammatitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gammatitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_ggsigma_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gammatitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')   ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gammatitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
		}
		}
		}
		
		else if "`distribution'" == "weibull" {
		
		constraint define 1 [p]_cons=1
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ggsigma_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`weibulltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`weibulltitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_ggsigma_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`weibulltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')   ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`weibulltitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
		}
		
		}
		}
		
		else if "`distribution'" == "ggamma" {
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ggsigma_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ggtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_ggsigma_group
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')   ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ggtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
		}
		
		}
		}
		
		* Fits a Normal Distribution if none specified
		else{
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_normal_group
			
			di " "
			di as txt "Fitting Full model with no constant:"
			
			
			if "`initial'" !=""{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial',copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			}
			
			else{
			
			local evaluator intllf_normal_group
			
			if "`initial'" !=""{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial', copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display			
			}
			}
		}		
		
}	
*Evaluates model if non-grouped data
else{	
		if "`distribution'" == "normal" {
		
			if "`noconstant'" != ""{
			
			local evaluator intllf_normal
			
			di " "
			di as txt "Fitting Full model with no constant:"
			
			
			if "`initial'" !=""{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial',copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			}
			
			else{
			local evaluator intllf_normal
			
			if "`initial'" !=""{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial', copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display			
			}
			}
		}
		
		else if "`distribution'" == "lnormal" {
		
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ln
			
			di " "
			di as txt "Fitting Full model with no constant:"
			
			if "`initial'" !=""{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial',copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			initial(`initial') `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			}
			
			else{
			local evaluator intllf_ln
			
			if "`initial'" !=""{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial', copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`lntitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display			
			}
			}
		}
		
		else if "`distribution'" == "t" {
		
			constraint define 1 [p]_cons=2
			constraint define 2 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"

			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
		
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "st" {
		
			constraint define 1 [p]_cons=2
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"

			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
		
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "gt" {
		
		constraint define 1 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "sgt" {
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			  [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			 (lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sgttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4)
			}
			
			}
			
			else{
			
			local evaluator intllf_sgt_condition
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q:`qvars') (lambda: `lambdavars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(4)
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sgttitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(4) 
			
			}
			

			}
		}
		
		else if "`distribution'" == "laplace" {
		
		constraint define 1 [p]_cons=1
		constraint define 2 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`laplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`laplacetitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`laplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints(1 2) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`laplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		
		else if "`distribution'" == "slaplace" {
		
		constraint define 1 [p]_cons=1
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`slaplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`slaplacetitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`slaplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`slaplacetitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		
		else if "`distribution'" == "ged" {
		
		constraint define 1 [lambda]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gedtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints (1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "sged" {
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`sgedtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`sgedtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "snormal" {
			
			constraint define 1 [p]_cons=2
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`snormaltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars')///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`snormaltitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display, neq(3) 
			}
			
			}
			
			else{
			
			local evaluator intllf_sged
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_normal (mu: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' `constraints'
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`snormaltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display, neq(3) 
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (lambda: `lambdavars') ///
			 [`weight'`exp'] if `touse' ==1 , maximize missing search(norescale) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`snormaltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat')

			ml display, neq(3) 
			
			}
			

			}
			
		}
		
		else if "`distribution'" == "br12" {
		
		constraint define 1 [p]_cons=1

			if "`noconstant'" != ""{
			
			local evaluator intllf_gb2exp
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(`initial' coeff2,copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br12title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			qui ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br12title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_gb2exp
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br12title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			
			di " "
			di as txt "Fitting Full model:"
						
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br12title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
				}
			}
			
			}
		
		else if "`distribution'" == "br3" {
		
		constraint define 1 [q]_cons=1

			if "`noconstant'" != ""{
			
			local evaluator intllf_gb2exp
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(`initial' coeff2,copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br3title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			qui ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br3title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_gb2exp
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`br3title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			
			di " "
			di as txt "Fitting Full model:"
						
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`br3title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
				}
			}
			
			}
			
		else if "`distribution'" == "gb2" {

			if "`noconstant'" != ""{
			
			local evaluator intllf_gb2exp
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')   ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(`initial' coeff2,copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			qui ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gb2title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_gb2exp
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			
			di " "
			di as txt "Fitting Full model:"
						
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars') (q: `qvars')  ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gb2title') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
				}
			}
			
			}
		
		else if "`distribution'" == "gamma" {

		constraint define 1 [lnsigma]_cons=0
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ggsigma
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gammatitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gammatitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_ggsigma
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`gammatitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')   ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`gammatitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
		}
		
		}
		}
		
		else if "`distribution'" == "weibull" {
		
		constraint define 1 [p]_cons=1
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ggsigma
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`weibulltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`weibulltitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_ggsigma
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`weibulltitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')   ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`weibulltitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
		}
		
		}
		}
		
		else if "`distribution'" == "ggamma" {
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_ggsigma
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs', noconstant)  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			[`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars' ) (p: `pvars') ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ggtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			
			}
			
			else{
			
			local evaluator intllf_ggsigma
			
			if "`initial'" != ""{
			
			*This portion here first evaluates the beta coefficients to get an estimate to pass it in later 
			* as start values
			
			qui ml model lf intllf_ln (delta: `depvar1' `depvar2'= `regs')  (lnsigma: `sigmavars' ) [`weight'`exp'] ///
			if `touse' ==1,missing search(norescale) maximize `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' `vce'  ///
			`robust' `cluster' 
			
			matrix coeff = e(b)
			
			if `nsigma' == 0{
				matrix coeff2 = coeff
				*matrix coeff2 = coeff[1..., 1..`nregs']
			}
			
			else{
				matrix coeff2 = coeff
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')  ///
			 [`weight'`exp'] if `touse' ==1 , maximize continue missing search(norescale) init(coeff2 `initial',copy) ///
			`constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  `nonrtolerance' /// 
			 `tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///  
			`robust' `cluster' `svy' repeat(`repeat') 
			
			ml display
			
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars' ) (p: `pvars')   ///
			[`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			 `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `tolerance' `ltolerance' `nrtolerance' `nonrtolerance'   ///  
			`showtolerance' title(`ggtitle') `vce'  /// 
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
		}
		
		}
		}
		
		* Fits a Normal Distribution if none specified
		else{
			
			if "`noconstant'" != ""{
			
			local evaluator intllf_normal
			
			di " "
			di as txt "Fitting Full model with no constant:"
			
			
			if "`initial'" !=""{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial',copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs',noconstant)  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance' `nonrtolerance' ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			}
			}
			
			else{
			
			local evaluator intllf_normal
			
			if "`initial'" !=""{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			init(`initial', copy) `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display
			
			}
			
			else{
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (lnsigma: ///
			`sigmavars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize /// 
			 `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'   ///
			`tolerance' `ltolerance' `nrtolerance' `nonrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			ml display			
			}
			}
		}		
		
}	
	*******************************************************************************************************************
		
		*Can't remember what the following section is for so I'm going to comment it out 
		
		/*
			if "`noconstant'" != ""{ 
			
			if "`distribution'" == "ln" | "`distribution'" == "lnormal" {
			
			local evaluator intllf_ln
			
			quietly ml model lf `evaluator' (mu: `depvar1' `depvar2' = )  /// 
			(lnsigma: ) [`weight'`exp'] if `touse' ==1 , missing search(on)  /// 
			maximize initial(`initial') `constraints' `technique'  `difficult'  ///
			`iterate' `log' `trace' `gradient' `showstep' `hessian'  ///
			`showtolerance' `tolerance' `ltolerance' `nrtolerance' title(`lntitle') ///
			`vce'  `robust' `cluster' `svy' repeat(`repeat')
			
			quietly ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  /// 
			(lnsigma: `sigmavars' ) [`weight'`exp'] if `touse' ==1 , missing search(on) continue /// 
			maximize initial(`initial') `constraints' `technique'  `difficult'  ///
			`iterate' `log' `trace' `gradient' `showstep' `hessian'  ///
			`showtolerance' `tolerance' `ltolerance' `nrtolerance' title(`lntitle') ///
			`vce'  `robust' `cluster' `svy' repeat(`repeat')
			
			}
			
			else if "`distribution'" == "gg" {
			
			local evaluator intllf_ggsigma
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = )   ///
			(lnsigma: ) (p: ) [`weight'`exp']  ///
			 if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(lnsigma: `sigmavars') (p: `pvars') [`weight'`exp']   ///
			 if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			}
			
			else if "`distribution'" =="sged"{
			
			local evaluator intllf_sged
			
			quietly ml model lf `evaluator' (mu: `depvar1' `depvar2' = )   ///
			(lambda: ) (lnsigma:  ) (p: ) [`weight'`exp']  ///
			 if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			quietly ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lambda: `lambdavars') (lnsigma: `sigmavars' ) (p: `pvars') [`weight'`exp']   ///
			 if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			}
			
			else if "`distribution'" == "gb2"{

			local evaluator intllf_gb2exp
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars') (p: `pvars') (q:  ///
			`qvars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = )   ///
			(lnsigma:  ) (p: ) (q:  ///
			) [`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			}
			
			else if "`distribution'" == "br3"{

			local evaluator intllf_gb2exp
			
			constraint define 1 [q]_cons=1
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars') (p: `pvars') (q:  ///
			`qvars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			initial(`initial') constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = )   ///
			(lnsigma:  ) (p: ) (q:  ///
			) [`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			}
			
			else if "`distribution'" == "br12"{

			local evaluator intllf_gb2exp
			
			constraint define 1 [q]_cons=1
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs', noconstant)   ///
			(lnsigma: `sigmavars') (p: `pvars') (q:  ///
			`qvars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize ///
			initial(`initial') constraints(1)`technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			quietly ml model lf `evaluator' (delta: `depvar1' `depvar2' = )   ///
			(lnsigma:  ) (p: ) (q:  ///
			) [`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') constraints(1) `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy' repeat(`repeat')
			
			}
			
			}
			
			*/
		
		mat betas = e(b) //coefficient matrix
		
		*Transforming lambda
		if "`distribution'" == "sgt" | "`distribution'" == "sged" | "`distribution'" == "gt" | ///
			"`distribution'" == "slaplace" | "`distribution'" == "laplace" | "`distribution'" == "st" | ///
			"`distribution'" == "ged" |"`distribution'" == "t" {
		
		di "{res:lambda}       {c |}"
		
		mat A = betas[1,"lambda:_cons"]
		mat A_SD = _se["lambda:_cons"]
		scalar alpha_c = A[1,1]
		scalar alpha_sd = A_SD[1,1]
		scalar lambda = (exp(alpha_c)- 1)/(exp(alpha_c)+ 1)
		scalar lambda_sd = alpha_sd*(2*exp(alpha_c)) / (exp(alpha_c) + 1)^2
		scalar zscore = lambda / lambda_sd
		scalar p = normal(-zscore)
		scalar llimit = lambda + 1.96*lambda_sd 
		scalar ulimit = lambda - 1.96*lambda_sd
		
		if lambda == 0{
		
		table_line_zero " _cons"  
		di as text "{hline 13}{c BT}{hline 64}"
		
		}
		
		else{
		table_line "_cons" lambda lambda_sd zscore p ulimit llimit 
		di as text "{hline 13}{c BT}{hline 64}"
		ereturn scalar lambda = lambda
		}
		
		}
		
		*Find the Conditional expected value at specified level
		if "`distribution'" == "gb2" | "`distribution'" == "ggamma" ///
		| "`distribution'" == "weibull" | "`distribution'" == "gamma" ///
		| "`distribution'" == "br12" | "`distribution'" == "br3" ///
		| "`distribution'" == "lnormal" ///
		| "`distribution'" == "st" | "`distribution'" == "gt" ///
		| "`distribution'" == "t" | "`distribution'" == "ged" ///
		| "`distribution'" == "sgt" | "`distribution'" == "sged" ///
		| "`distribution'" == "slaplace" | "`distribution'" == "laplace" ///
		| "`distribution'" == "snormal" {
			
			mat mid_Xs = 1
			if "`eyx'" == ""{
				local eyx "mean"
				di "{res:`eyx'}         {c |}"
			}
			else if "`eyx'" == "mean" {
			di "{res:`eyx'}          {c |}" 
			}
			else if "`eyx'" == "p50" | "`eyx'" == "p10" | "`eyx'" == "p25" | ///
			        "`eyx'" == "p75" | "`eyx'" == "p90" | "`eyx'" == "p95" | ///
					"`eyx'" == "p99" | "`eyx'" == "min" | "`eyx'" == "max" {		
			
					di "{res:`eyx'}          {c |}"
				}
			else if "`eyx'" == "p1" | "`eyx'" == "p5" {
				di "{res:`eyx'}             {c |}"
				}
			else{
				di as err "Not a valid option for eyx"
				exit 498
				}
			
			
			quietly foreach x in `regs' {
				sum `x', detail
				scalar mid_ = r(`eyx')
				mat mid_Xs = mid_Xs, mid_
			}
			mat sigma = betas[1,"lnsigma:_cons"]
			scalar sigma = sigma[1,1]
			
			if "`distribution'" == "gb2" | "`distribution'" == "br12" | "`distribution'" == "br3" {
			
				mat deltas = betas[1,"delta:"]
				mat deltas = deltas'
				mata: st_matrix("deltas", flipud(st_matrix("deltas"))) //flips matrix around
		                        									// to conform with Xs									
				mat p = betas[1,"p:_cons"]
				scalar p = p[1,1]
				mat q = betas[1,"q:_cons"]
				scalar q = q[1,1]
				mat xbeta = mid_Xs*deltas
				scalar xbeta = xbeta[1,1]
				mat expected = exp(xbeta)*( (exp(lngamma(p+exp(sigma)))*exp(lngamma(q-exp(sigma))))/  ///
											( exp(lngamma(p))*exp(lngamma(q))))
			}
			
			
			if "`distribution'" == "sgt" | "`distribution'" == "st" | ///
			"`distribution'" == "gt" | "`distribution'" == "t" {
			
				mat mu = betas[1,"mu:"]
				mat mu = mu'
				mata: st_matrix("mu", flipud(st_matrix("mu"))) //flips matrix around
																		// to conform with Xs										
				mat xbeta = mid_Xs*mu     																
				mat p = betas[1,"p:_cons"]
				scalar p = p[1,1]
				mat q = betas[1,"q:_cons"]
				scalar q = q[1,1]
				mat A = betas[1,"lambda:_cons"]
				scalar alpha_c = A[1,1]
				scalar lambda = (exp(alpha_c)- 1)/(exp(alpha_c)+ 1)
				mat sigma = betas[1,"lnsigma:_cons"]
				scalar sigma = sigma[1,1]
				scalar xbeta = xbeta[1,1]
				mat expected = xbeta + 2*lambda*exp(sigma)*((q^(1/p))*(exp(lngamma(2/p) ///
				+lngamma(q-(1/p)) - lngamma((1/p)+q))/exp(lngamma(1/p) ///
				+lngamma(q)) - lngamma((1/p)+q)) ) 
			}
			
			if "`distribution'" == "sged" | "`distribution'" == "ged" | ///
				"`distribution'" == "slaplace" | "`distribution'" == "laplace" | ///
				"`distribution'" == "snormal" {
			
				mat mu = betas[1,"mu:"]
				mat mu = mu'
				mata: st_matrix("mu", flipud(st_matrix("mu"))) //flips matrix around
																		// to conform with Xs										
				mat xbeta = mid_Xs*mu     																
				mat p = betas[1,"p:_cons"]
				scalar p = p[1,1]
				mat A = betas[1,"lambda:_cons"]
				scalar alpha_c = A[1,1]
				scalar lambda = (exp(alpha_c)- 1)/(exp(alpha_c)+ 1)
				mat sigma = betas[1,"lnsigma:_cons"]
				scalar sigma = sigma[1,1]
				scalar xbeta = xbeta[1,1]
				mat expected = xbeta + 2*lambda*exp(sigma)*(exp(lngamma(2/p)) /exp(lngamma(1/p)))		
			}
			
			if "`distribution'" == "ggamma" | "`distribution'" == "gamma" | "`distribution'" == "weibull" {
			
				mat deltas = betas[1,"delta:"]
				mat deltas = deltas'
				mata: st_matrix("deltas", flipud(st_matrix("deltas"))) //flips matrix around
																		// to conform with Xs										
				mat p = betas[1,"p:_cons"]
				scalar p = p[1,1]
				mat xbeta = mid_Xs*deltas
				scalar xbeta = xbeta[1,1]
				mat expected = exp(xbeta)*( (exp(lngamma(p+exp(sigma))))/  ///
											( exp(lngamma(p))))
			}
			
			if "`distribution'" == "lnormal" {
				
				mat mu = betas[1,"mu:"]
				mat mu = mu'
				mata: st_matrix("mu", flipud(st_matrix("mu"))) //flips matrix around
																		// to conform with Xs										
				mat xbeta = mid_Xs*mu
				scalar xbeta = xbeta[1,1]
				mat expected = exp(xbeta + (exp(sigma)^2/2))
			}
					
			scalar eyx = expected[1,1]
			table_line "E[Y|X]" eyx 
			di as text "{hline 13}{c BT}{hline 64}"
			ereturn scalar eyx = eyx
		}
		
		if "`distribution'" == "normal" | "`distribution'" == "" {
			mat sigma = betas[1,"lnsigma:_cons"]
			scalar sigma = sigma[1,1]
			}

		qui count
		scalar numobs = r(N)
		scalar numvars = e(df_m)
		scalar logll = e(ll)
	
		if "`noconstant'" != ""{
			if "`distribution'" == "gb2"{
				if p == 1{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
				else if q == 1{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
				else{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
			}
			else if "`distribution'" == "lnormal" {
				scalar BIC = ln(numobs)*(numvars + 1)-2*logll
				scalar AIC = 2* (numvars + 1) - 2* logll
			}
			else if "`distribution'" == "ggamma" {

				if p == 1{
				scalar BIC = ln(numobs)*(numvars + 1)-2*logll
				scalar AIC = 2* (numvars + 1) - 2* logll
				}
				else if sigma == 0{
				scalar BIC = ln(numobs)*(numvars + 1)-2*logll
				scalar AIC = 2* (numvars + 1) - 2* logll
				}
				else{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
			}
			else if "`distribution'" == "sgt" {
				if p == 2{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
				else if lambda == 0{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars +3) - 2* logll
				}
				else if lambda == 0 & p ==2{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
				else{
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars + 4) - 2* logll
				}
			}
			else if "`distribution'" == "sged" {
				if lambda == 0{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
				else {
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars +3) - 2* logll
				}
			}
			
			else if "`distribution'" == "st" {
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars +3) - 2* logll
			}
			
			else if "`distribution'" == "gt" {
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars +3) - 2* logll
			}
			
			else if "`distribution'" == "ged" {
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
			}
			
			else if "`distribution'" == "weibull" {
				scalar BIC = ln(numobs)*(numvars + 1)-2*logll
				scalar AIC = 2* (numvars + 1) - 2* logll
			}
			
			else if "`distribution'" == "gamma" {
				scalar BIC = ln(numobs)*(numvars + 1)-2*logll
				scalar AIC = 2* (numvars + 1) - 2* logll
			}
			
			else if "`distribution'" == "t" {
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
			}
			
			else{
				scalar BIC = ln(numobs)*(numvars + 1)-2*logll
				scalar AIC = 2* (numvars + 1) - 2* logll
			}
		}
		else
		{
			if "`distribution'" == "gb2" {
				if p == 1{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
				else if q == 1{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
				else{
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars + 4) - 2* logll
				}
				
			}
			else if "`distribution'" == "lnormal" {
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
			}
			else if "`distribution'" == "ggamma" {

				if p == 1{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
				else if sigma == 0{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
				}
				else{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
			}
			else if "`distribution'" == "sgt" {
				if p == 2{
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars + 4) - 2* logll
				}
				else if lambda == 0{
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars + 4) - 2* logll
				}
				else if lambda == 0 & p ==2{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
				else{
				scalar BIC = ln(numobs)*(numvars + 5)-2*logll
				scalar AIC = 2* (numvars + 5) - 2* logll
				}
			}
			else if "`distribution'" == "sged" {
				if lambda == 0{
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
				}
				else {
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars + 4) - 2* logll
				}
			}
			
			else if "`distribution'" == "st" {
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars +4) - 2* logll
			}
			
			else if "`distribution'" == "gt" {
				scalar BIC = ln(numobs)*(numvars + 4)-2*logll
				scalar AIC = 2* (numvars +4) - 2* logll
			}
			
			else if "`distribution'" == "ged" {
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
			}
			
			else if "`distribution'" == "weibull" {
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
			}
			
			else if "`distribution'" == "gamma" {
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
			}
			else if "`distribution'" == "t" {
				scalar BIC = ln(numobs)*(numvars + 3)-2*logll
				scalar AIC = 2* (numvars + 3) - 2* logll
			}
			
			
			else{
				scalar BIC = ln(numobs)*(numvars + 2)-2*logll
				scalar AIC = 2* (numvars + 2) - 2* logll
			}
		ereturn scalar bic = BIC
		ereturn scalar aic = AIC
		}
				
		* Displays alternate notation for GB2, SGT trees
		if "`distribution'" == "gb2" | "`distribution'" == "br3" | "`distribution'" == "br12" | ///
			"`distribution'" == "gamma" | "`distribution'" == "ggamma" | "`distribution'" == "weibull" {
			di as txt "Alternate notation: "
			di as txt "a: " 1/exp(sigma)
			di as txt "beta: " exp(xbeta)
		}
		
		else if "`distribution'" == "lnormal" | "`distribution'" == "sgt" | ///
			"`distribution'" == "sged" | "`distribution'" == "slaplace" | "`distribution'" == "laplace" | ///
			"`distribution'" == "snormal" | "`distribution'" == "st" | "`distribution'" == "ged" | ///
			"`distribution'" == "gt" | "`distribution'" == "st" | "`distribution'" == "t" | ///
			"`distribution'" == "normal" {
			di as txt "Alternate notation: "
			di as txt "sigma: " exp(sigma)
		}
		
		else {
			di as txt "Alternate notation: "
			di as txt "sigma: " exp(sigma)
		}
		
		* Displays BIC and AIC values
		if "`bicaic'" != "" {
		    * only for non-grouped data
		    if "`frequency'" == "" {	
				di _n(0)
				di as txt "BIC: " BIC
				di as txt "AIC: " AIC
			}
			else {
			    di _n(0)
			    di as txt "Option BICAIC not compatible with grouped data"
			}
		}
		
		* Reports gini coefficient
		if "`gini'" != "" {
			di _n(0)
			if "`distribution'" == "weibull" {
				di as txt "Gini coefficient: " 1-(.5^(exp(sigma)))
			}
			if "`distribution'" == "gamma" {
				di as txt "Gini coefficient: " exp(lngamma(p + .5)) / (exp(lngamma(p + 1)) * sqrt(_pi))
			}			
			else if "`distribution'" == "br3" {
			    di as txt "Gini coefficient: " [exp(lngamma(p)) * exp(lngamma(2*p + exp(sigma)))] ///
												/ [exp(lngamma(p + exp(sigma))) * exp(lngamma(2*p))] - 1
			}
			else if "`distribution'" == "br12" {
			    di as txt "Gini coefficient: " 1 - [exp(lngamma(q)) * exp(lngamma(2*q - exp(sigma)))] ///
												/ [exp(lngamma(q - exp(sigma))) * exp(lngamma(2*q))]
			}
		}
		
		* plots pdf
		if "`plot'" != "" {
		    local range = "`plot'"
			local rangelow = min(`range')
			local rangehigh = max(`range')

		    if "`distribution'" == "gb2" | "`distribution'" == "br3" | "`distribution'" == "br12" {
				
				local a = 1/exp(sigma)
				local b = exp(xbeta)
				local p = p
				local q = q
				local B =  exp(lngamma(`p')+lngamma(`q')-lngamma(`p' + `q'))

				twoway function ///
					y = [abs(`a') * (x / `b')^(`a' * `p')] / [x * (`B') * (1 + (x / `b')^`a')^(`p' * `q')], ///
					range(`rangelow' `rangehigh')
			}
			
			else if "`distribution'" == "ggamma" | "`distribution'" == "gamma" | "`distribution'" == "weibull" {
			    
			    local a = 1/exp(sigma)
				local b = exp(xbeta)
				local p = p
				local G = exp(lngamma(`p'))
				
				***
				
				di `a'
				di `b'
				di `p'
				di `q'
				***
				
				twoway function ///
					y = (abs(`a')*(x/`b')^(`a'*`p')*exp(-(x/`b')^`a'))/(x*`G'), ///
					range(`rangelow' `rangehigh')
			}
			
			else if "`distribution'" == "sgt" | "`distribution'" == "st" | ///
					"`distribution'" == "t" | "`distribution'" == "gt"  {
				
				local lambda = lambda
				local s = exp(sigma)
				local p = p
				local q = q
				local B =  exp(lngamma(1/`p')+lngamma(`q')-lngamma(1/`p'+`q'))

				twoway function ///
					y = (`p')/[(2*`s'*`q'^(1/`p')*`B')*(1 + (abs(x)^`p')/(`q'*`s'^`p'*(1 + `lambda'*sign(x))^`p'))^(`q' + 1/`p')], ///
					range(`rangelow' `rangehigh')
			}
			
			else if "`distribution'" == "sged" | "`distribution'" == "ged" | ///
					"`distribution'" == "slaplace" | "`distribution'" == "laplace" | ///
					"`distribution'" == "snormal" {

				local s = exp(sigma)
				local lambda = lambda
				local p = p
				local G = exp(lngamma(1/`p'))

				twoway function ///
					y = [`p' * exp(-(x^`p' / ((1 + `lambda' * sign(x))^`p' * `s'^`p')))] / [2 * `s' * `G'], ///
					range(`rangelow' `rangehigh')
			}
						
			else if "`distribution'" == "lnormal" {
				
				local mu = mu[1,1]				
				local sigma = exp(sigma)

				twoway function ///
				y = [exp(-(ln(x) - `mu')^2 / 2 * `sigma'^2)] / [sqrt(2 * c(pi)) * x * `sigma'], ///
				range(`rangelow' `rangehigh')
			}
			
			else if "`distribution'" == "normal" | "`distribution'" == "" {
			    
				mat mu = betas[1,"mu:_cons"]
				mat sigma = betas[1,"lnsigma:_cons"]
				
				local mu = mu[1,1]
				local sigma = exp(sigma[1,1])
				
				twoway function ///
					y = normalden(x, `mu', `sigma'), ///
					range(`rangelow' `rangehigh')
			}
		}
		
		
		*Observation type count for interval regression
		if "`frequency'" == ""{
			noi di " "
			if `nleft' != 1{
				noi di as txt " {res:`nleft'} left-censored observations" 
			}
			if `nleft' == 1{
				noi di as txt " {res:`nleft'} left-censored observation" 
			}
			if `nuncensored' != 1{
				noi di as txt " {res: `nuncensored'} uncensored observations" 
			}
			if `nuncensored' == 1{
				noi di as txt " {res:`nuncensored'} uncensored observation" 
			}
			if `nright' != 1{
				noi di as txt " {res:`nright'} right-censored observations" 
			}
			if `nright' == 1{
				noi di as txt " {res:`nright'} right-censored observation" 
			}
			if `ninterval' != 1{
				noi di as txt " {res:`ninterval'} interval observations" 
			}
			if `ninterval' == 1{
				noi di as txt " {res:`ninterval'} interval observation" 
			}
		}
		*Observation type count for grouped regression
		if "`frequency'" ~= ""{
			
			noi di " "
			noi di as txt " {res: `total'} groups"
			if `nleft' != 1{
				noi di as txt " {res:`nleft'} left-censored groups" 
			}
			if `nleft' == 1{
				noi di as txt " {res:`nleft'} left-censored group" 
			}
			if `nuncensored' != 1{
				noi di as txt " {res: `nuncensored'} uncensored groups" 
			}
			if `nuncensored' == 1{
				noi di as txt " {res:`nuncensored'} uncensored group" 
			}
			if `nright' != 1{
				noi di as txt " {res:`nright'} right-censored groups" 
			}
			if `nright' == 1{
				noi di as txt " {res:`nright'} right-censored group" 
			}
			if `ninterval' != 1{
				noi di as txt " {res:`ninterval'} interval groups" 
			}
			if `ninterval' == 1{
				noi di as txt " {res:`ninterval'} interval group" 
			}
		}
		qui ereturn list
		}
	
end

*program drop table_line
capture program drop table_line
program table_line
	args vname coef se z p 95l 95h
	if (c(linesize) >= 100){
		local abname = "`vname'"
		}
	else if (c(linesize) > 80){
	local abname = abbrev("`vname'", 12+(c(linesize)-80))
	}
	else{
	local abname = abbrev("`vname'", 12)
	}
	local abname = abbrev("`vname'",12)
	display as text %12s "`abname'" " { c |}" /*
	*/ as result /*
	*/ "    " %8.0g `coef' "  " /*
	*/ %9.0g `se' "   " %03.2f `z' "   " /*
	*/ %04.3f `p' "    " %9.0g `95l' "    " /*
	*/ %9.0g `95h'
end

capture program drop table_line_zero
program table_line_zero
	args vname 
	if (c(linesize) >= 100){
		local abname = "`vname'"
		}
	else if (c(linesize) > 80){
	local abname = abbrev("`vname'", 12+(c(linesize)-80))
	}
	else{
	local abname = abbrev("`vname'", 12)
	}
	local abname = abbrev("`vname'",12)
	display as text %12s "`abname'" " { c |}" /*
	*/ as result /*
	*/ "          " "0  {sf:(constrained)}" " " 
end

capture program drop ParseHet
program ParseHet, rclass
        syntax varlist(fv ts numeric) [, noCONStant]
        return local varlist "`varlist'"
        return local constant `constant'
end

version 13.0
mata:
matrix function flipud(matrix X)
{
return(rows(X)>1 ? X[rows(X)..1,.] : X)
}
end
