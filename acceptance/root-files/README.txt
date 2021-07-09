## Acceptance is generally calculated for testing with cuts:
## 	sim:	PhiPQ!=-9999, TargType!=-9999, pid==211; mc_PhiPQ!=-9999, mc_TargType==2, mc_pid==211.
## When adding Nphe cuts, these are implemented only for reconstructed data (i.e. not mc_<> variables).
## Easy cuts (S.Moran): Q2, Xb, Zh, Pt2, PhiPQ as seen in table 4.13 (p. 80) of his thesis

Acceptance_<target>1_nNphe.root			(4)	: Acceptance calculated without Nphe cut. (and simple cuts)
Acceptance_<target>1_<l/g>25.root		(5)	: Acceptance calculated with Nphe cut, lower or greater than 25.
Acc_Fe1_CTest_l25.root				(1)	: First attempt to make use of Closure Test (CT) with Acceptance Nphe<25
Acc_<target>1_cuts.root				(4)	: Acceptance with easy cuts and Nphe(hadrons)>25
Acc_Fe1_cuts2.root				(1)	: Acceptance with easy cuts and Nphe(hadrons)<25
