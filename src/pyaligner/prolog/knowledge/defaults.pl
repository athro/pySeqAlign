:- assert(gapDefault(-1.0)).
:- assert(gapChar('$gap')).
:- assert(gapChar('real_gap')).
:- assert(learningRate(0,1.0)).

assign(X,V) :-
	Old =..[X,_], retract(Old),
	New =..[X,V], assert(New).

dist(sym,atomDistance,nc,0,0,_,Dist):-	gapDefault(Dist).
dist(sym,atomDistance,nc,0,_,0,Dist):-	gapDefault(Dist).

:- consult('aminoAcids.pl'),consult('blossum_50.pl').


