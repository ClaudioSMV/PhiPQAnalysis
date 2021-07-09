## Acceptance is generally calculated for testing with cuts:
## 	sim:	PhiPQ!=-9999, TargType!=-9999, pid==211; mc_PhiPQ!=-9999, mc_TargType==2, mc_pid==211.
## When adding Nphe cuts, these are implemented only for reconstructed data (i.e. not mc_<> variables).
## Easy cuts (S.Moran): Q2, Xb, Zh, Pt2, PhiPQ as seen in table 4.13 (p. 80) of his thesis.
##

NAME						N°FILES	: DESCRIPTION

Nphevs_1<file>.root				(2)	: TH2F Nphe vs (Q2, Nu, Zh, Pt2, PhiPQ) (pid, TargType and Nphe [0,200] cuts, Fe)
Nphevs_simN<2/5>.root				(2)	: Same as above, but Xb instead of Nu and cuts of S.Moran + <2/5> binning Nphe, Fe
Nphevs_<target><file>.root			(6)	: Same as above, 2 Nphe binning, corrected binning for other variables. <sim/data>
