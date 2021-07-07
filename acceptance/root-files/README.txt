## Acceptance is generally calculated for testing with cuts:
## 	sim:	PhiPQ!=-9999, TargType!=-9999, pid==211; mc_PhiPQ!=-9999, mc_TargType==2, mc_pid==211.
## When adding Nphe cuts, these are implemented only for reconstructed data (i.e. not mc_<> variables).
##
## Analysis is generally calculated for testing with cuts:
##	data:	TargType correct, pid==211, (Nphe>0), Nu>=2, Q2>1.35, (Q2<1.82)
##

Analysis_OLD.root				(1)	: Old analysis, not clear.
Analysis_Fe1_old<n>.root + Analysis_Fe1_2.root	(3)	: Old analysis files (check what differences do they have).
Analysis_Fe1_<l/g>25.root			(2)	: Correction of data using <l/g> acceptance (above) and same cut in Nphe.

Acceptance_<target>1_nNphe.root			(4)	: Acceptance calculated without Nphe cut. (and simple cuts)
Acceptance_<target>1_<l/g>25.root		(5)	: Acceptance calculated with Nphe cut, lower or greater than 25.

Nphevs_1<file>.root				(2)	: TH2F Nphe vs (Q2, Nu, Zh, Pt2, PhiPQ) (pid, TargType and Nphe [0,200] cuts, Fe)
Nphevs_simN<2/5>.root				(2)	: Same as above, but Xb instead of Nu and cuts of S.Moran + <2/5> binning Nphe, Fe
Nphevs_<target><file>.root			(6)	: Same as above, 2 Nphe binning, correct binning for other variables. <sim/data>

Acc_Fe1_CTest_l25.root				(1)	: First attempt to make use of Closure Test (CT) with Acceptance Nphe<25
