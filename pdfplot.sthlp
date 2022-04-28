{smcl}
{* *! version 1.0  28apr2022}{...}
{viewerjumpto "Syntax" "examplehelpfile##syntax"}{...}
{viewerjumpto "Description" "examplehelpfile##description"}{...}
{viewerjumpto "Examples" "examplehelpfile##examples"}{...}
{title:Title}

{phang}
{bf:pdfplot} {hline 2} Plot a Distribution's Probability Density Function


{marker syntax}{...}
{title:Syntax}

{p 8 20 2}
{cmdab:pdfplot}
{it:dist_name}
{cmd:,}
{it: params({help numlist:numlist})}
{it: {help twoway_function:range()}}


{pstd}
Supported distributions come from the Skewed Generalized t and Generlized Beta of the Second Kind families of distributions, namely: 

		Distribution{space 25}{it:dist_name}{space 6}{it:params}
		{hline 71}
		Normal{space 31}normal{space 9}mu, sigma
		Lognormal{space 28}lnormal{space 8}mu, sigma
		Skewed Normal{space 24}snormal{space 8}sigma, lambda
		Laplace{space 30}laplace{space 8}sigma
		Skewed Laplace{space 23}slaplace{space 7}sigma, lambda
		Generalized Error{space 20}ged{space 12}sigma, p
		Skewed Generalized Error{space 13}sged{space 11}sigma, p, lambda
		t{space 36}t{space 14}sigma, q
		Generalized t{space 24}gt{space 13}sigma, p, q
		Skewed t{space 29}st{space 13}sigma, q, lambda
		Skewed Generalized t{space 17}sgt{space 12}sigma, p, q, lambda
		Weibull{space 30}weibull{space 8}a, b
		Gamma{space 32}gamma{space 10}b, p
		Generalized Gamma{space 20}ggamma{space 9}a, b, p
		Burr 3{space 31}br3{space 12}a, b, p
		Burr 12{space 30}br12{space 11}a, b, q
		Generalized Beta of the Second Kind{space 2}gb2{space 12}a, b, p, q
		{hline 71}

		
{marker description}{...}
{title:Description}

{pstd}
{cmd:pdfplot} is a tool to efficiently plot the PDF of one of the above distributions with chosen parameter values. It uses {help twoway_function} to do so. {cmd:pdfplot} was originally developed as an auxialliary function to {cmd:{help gintreg}}, which performs maxiumum likelihood estimation on any of the same set of accepted distributions and is compatible with numerous forms of data, including interval, truncated, and frequency in addition to the typical form.


{marker examples}{...}
{title:Examples}

{pstd}Plotting a Standard Normal PDF{p_end}
{phang2}{cmd:. pdfplot normal, params(0 1) range(-4 4)}{p_end}

{pstd}Plotting a Generalized Gamma PDF with arbitrary parameters{p_end}
{phang2}{cmd:. pdfplot ggamma, params(10 20 .3) range(0 30)}{p_end}

{pstd}Plotting a Skewed Generalized t PDF with an adjusted lambda after estimating with {cmd:gintreg}{p_end}
{phang2}{cmd:. webuse intregxmpl}{p_end}
{phang2}{cmd:. gintreg wage1 wage2, dist(sgt)}{p_end}
{phang2}{cmd:. pdfplot sgt, params(7.423 1.743 1.846 1) range(0 50)}{p_end}


{marker author}{...}
{title:Author}

{phang}
Authored by Jacob Triplett at Brigham Young University. For support contact Jacob at jtriplett100@gmail.com.
{p_end}