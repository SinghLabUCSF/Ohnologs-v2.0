## The core scripts to identify ohnologs


* Needs perl v-5.0 or above with these modules installed:     
	* Math::Combinatorics
	* Math::BigInt
	* Math::BigFloat
	
* The paths to the [parameter files](../2_ParameterFiles/) are relative path starting from this directory. Sample parameter files for [hsapiens](../2_ParameterFiles/hsapiens) 

> Sample commands: 
> 
>   `perl OHNOLOGS_SelfComparison_CL.pl '../2_ParameterFiles/hsapiens/self' '../3_SyntenyOutputFiles/hsapiens/self' 2>&1 | tee hsapiens_self.log`
>   
>   `perl OHNOLOGS_OutgroupComparison_CL.pl '../2_ParameterFiles/hsapiens/outgroup' '../3_SyntenyOutputFiles/hsapiens/outgroup' 2>&1 | tee hsapiens_og.log`