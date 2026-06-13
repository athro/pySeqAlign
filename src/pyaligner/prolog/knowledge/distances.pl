%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       %
%      distances for logical atoms      %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The predicat dist/6 is meant to be called witch instanciated atoms.
%
% The form is:
% dist(TypeOfDistance,NameOfDistance,Iteration,Atom1,Atom2,Distance)
%
%
%

%:- use(module(library(lists))).
:- source.

%:- dynamic gapDefault/1.
%:- dynamic gapChar/1.
:- dynamic x/1 .	      % this may be required in some Prologs  
:- assert(gapDefault(-1.0)).
:- assert(gapChar('$gap')).
:- assert(gapChar('real_gap')).
:- assert(learningRate(0,1.0)).

x(0).			% An initial value is required in this example

assign(X,V) :-
	Old =..[X,_], retract(Old),
	New =..[X,V], assert(New).

% Nienhuys-Cheng Distance


dist(dist,TypeOfDistance,NameOfDistance,0,AtomID1,AtomID2,Distance) :-
	example(AtomID1,Atom1),
	gapChar(Atom1),
	gapDefault(Distance).
dist(dist,TypeOfDistance,NameOfDistance,0,AtomID1,AtomID2,Distance) :-
	example(AtomID2,Atom2),
	gapChar(Atom2),
	gapDefault(Distance).
dist(dist,TypeOfDistance,NameOfDistance,0,AtomID1,AtomID2,Distance) :-
	example(AtomID1,Atom1),
	example(AtomID2,Atom2),
	distSub(TypeOfDistance,NameOfDistance,0,Atom1,Atom2,Distance).

dist(sym,TypeOfDistance,NameOfDistance,0,AtomID1,AtomID2,Distance) :-
	example(AtomID1,Atom1),
	gapChar(Atom1),
	gapDefault(Distance).
dist(sym,TypeOfDistance,NameOfDistance,0,AtomID1,AtomID2,Distance) :-
	example(AtomID2,Atom2),
	gapChar(Atom2),
	gapDefault(Distance).
dist(sym,TypeOfDistance,NameOfDistance,0,AtomID1,AtomID2,Distance) :-
	example(AtomID1,Atom1),
	example(AtomID2,Atom2),
	distSub(TypeOfDistance,NameOfDistance,0,Atom1,Atom2,DistanceI),
	Distance is 1.0-DistanceI.
	

distSub(atomDistance,nc,0,Atom,Atom,Dist):-
	!,Dist is 0.0.
%distSub(atomDistance,nc,0,A,_,Dist):-
%	gapChar(A),
%	gapDefault(Dist).
%distSub(atomDistance,nc,0,_,B,Dist):-
%	gapChar(B),
%	gapDefault(Dist).
distSub(atomDistance,nc,0,A,B,Dist) :-
	A =.. [PredA|AL],
	B =.. [PredB|BL],
	PredA == PredB,
	length(AL,Length),
	length(BL,Length),!,
	distSub_helper(atomDistance,nc,0,AL,BL,SumDist),
	Dist is 1.0/(2*Length)*SumDist.
distSub(atomDistance,nc,0,A,B,Dist) :-
	Dist is 1.0.

distSub_helper(atomDistance,nc,0,[],[],Dist) :-
	!,Dist is 0.0.
distSub_helper(atomDistance,nc,0,[A1|R1],[A2|R2],Dists) :-
	distSub(atomDistance,nc,0,A1,A2,DistHere),
	distSub_helper(atomDistance,nc,0,R1,R2,DistsThere),!,
	Dists is DistHere+DistsThere.


%distSub(atomDistance,nc,Iteration,Atom1,Atom2,Dist):-
%	IterBefore is Iteration-1,
%	distSub(atomDistance,nc,IterBefore,Atom1,Atom2,DistBefore),


% this delta/5 works on the prolog programs resulting from a tilde run.
delta(_, 0, _,_, 0.0).
delta(Num, Iteration, Atom1,Atom2, Delta) :-
	exampleC(Atom1,tag(What1,Tag1)),
	assert(word(Num,What1)),
	assert(tag(Num,What1,Tag1)),
	assert(q(Num,What1)),
	exampleC(Atom2,tag(What2,Tag2)),
	% only assert if different
	(What2 \== What1 ->
	    assert(word(Num,What1));
	    true),
	% only assert if different
	((Tag1 \== Tag2, What2 \== What1) ->
	    assert(tag(Num,What2,Tag2));
	    true),
	assert(t(Num,What2)),!,
	delta(Iteration, Num, [Delta]),
	retract(q(Num,What1)),
	retract(t(Num,What2)),
        (retract(word(Num,What1));true),
        (retract(word(Num,What2));true),
        (retract(tag(Num,What1,Tag1));true),
        (retract(tag(Num,What2,Tag2));true).

