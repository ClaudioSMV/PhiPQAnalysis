## Acceptance is generally calculated for testing with cuts:
## 	sim:	PhiPQ!=-9999, TargType!=-9999, pid==211; mc_PhiPQ!=-9999, mc_TargType==2, mc_pid==211.
## When adding Nphe cuts, these are implemented only for reconstructed data (i.e. not mc_<> variables).
## Easy cuts (S.Moran): Q2, Xb, Zh, Pt2, PhiPQ as seen in table 4.13 (p. 80) of his thesis

NAME						N°FILES	: DESCRIPTION
						(24)

ClosTestCSMV_<target>1.root			(4)	: Closure Test using all SMoran's cuts, Claudio's way, with correct Errors
ClosTestSM_<target>1.root			(4)	: Closure Test using all SMoran's cuts, SMoran's way, with correct Errors
ClosTestCSMV_<target>1_EX_<var>.root		(8)	: CT with all SMoran's cuts + bins in Q2&Nu, Claudio's way. <var>={Zh,Pt2}
ClosTestSM_<target>1_EX_<var>.root		(8)	: CT with all SMoran's cuts + bins in Q2&Nu, SMoran's way. <var>={Zh,Pt2}
