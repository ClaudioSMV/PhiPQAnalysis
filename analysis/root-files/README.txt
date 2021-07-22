## Analysis root files generated
#################################
## Analysis is generally calculated for testing with cuts:
##	data:	TargType correct, pid==211, (Nphe>0), Nu>=2, Q2>1.35, (Q2<1.82)
## Easy cuts (S.Moran): Q2, Xb, Zh, Pt2, PhiPQ as seen in table 4.13 (p. 80) of his thesis
##

NAME						N°FILES	: DESCRIPTION

Analysis_OLD.root				(1)	: Old analysis, not clear
Analysis_Fe1_old<n>.root + Analysis_Fe1_2.root	(3)	: Old analysis files (check what differences do they have)
Analysis_Fe1_<l/g>25.root			(2)	: Correction of data using <l/g> acceptance and same cut in Nphe
AnalysisSM_Fe1.root				(1)	: Correction using SMorán's implementation of both, acceptance and correction
AnalysisSM_Fe1_light.root			(1)	: Same as above but without all corrected bins saved in the file (lighter)
AnalysisCSMV_Fe1.root				(1)	: Correction using Claudio's implementation of both, acceptance and correction
