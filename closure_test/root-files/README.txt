## Acceptance is generally calculated for testing with cuts:
## 	sim:	PhiPQ!=-9999, TargType!=-9999, pid==211; mc_PhiPQ!=-9999, mc_TargType==2, mc_pid==211.
## When adding Nphe cuts, these are implemented only for reconstructed data (i.e. not mc_<> variables).
## Easy cuts (S.Moran): Q2, Xb, Zh, Pt2, PhiPQ as seen in table 4.13 (p. 80) of his thesis

NAME						N°FILES	: DESCRIPTION
						(4)

CT_PhiPQ_1 (folder)				(1)	: Closure Test using all SMoran's cuts, with correct Errors
CT_EX_1 (folder)				(1)	: CT with all SMoran's cuts + bins in Q2&Nu, Claudio's way. <var>={Zh,Pt2}
CT_PhiPQ_2 (folder)				(1)	: Same as EX_1, but with PhiPQ variable (maybe not well binned)
CT_EX_2 (folder)				(1)	: Zh, Pt2, PhiPQ with (now well understand) acceptance and corr in diff bins
