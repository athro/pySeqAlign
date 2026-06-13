example(0,'-').
example(1,'a').
example(2,'r').
example(3,'n').
example(4,'d').
example(5,'c').
example(6,'q').
example(7,'e').
example(8,'g').
example(9,'h').
example(10,'i').
example(11,'l').
example(12,'k').
example(13,'m').
example(14,'f').
example(15,'p').
example(16,'s').
example(17,'t').
example(18,'w').
example(19,'y').
example(20,'v').

getExampleIDs([],[]).
getExampleIDs([Haa|Taa],[Hid|Tid]) :-
	example(Hid,Haa),
	getExampleIDs(Taa,Tid).


%MAPFQSNKDL
%A = [m,a,p,f,q,s,n,k,d,l],
%S1 = [13,1,15,14,6,16,3,12,4,11] ?
%A = ['m', 'a', 'p', 'f', 'q', 's', 'n', 'k', 'd', 'l'], getExampleIDs(A,S1).
sequence(1,[13,1,15,14,6,16,3,12,4,11]).

% MLAPFEKTAAARSII
%getExampleIDs(['m', 'l', 'a', 'p', 'f', 'e', 'k', 't', 'a', 'a', 'a', 'r', 's', 'i', 'i'],Seq2).
sequence(2,[13,11,1,15,14,7,12,17,1,1,1,2,16,10,10]).

%getExampleIDs(['h', 'e', 'a', 'g', 'a', 'w', 'g', 'h', 'e', 'e'],Seq3)
sequence(3,[9,7,1,8,1,18,8,9,7,7]).
%getExampleIDs(['p', 'a', 'w', 'h', 'e', 'a', 'e'],Seq4)
sequence(4,[15,1,18,9,7,1,7]).

ppAAList([]) :- nl.
ppAAList([H|T]) :-
	format("~w",[H]),
	ppAAList(T).

convertToString(TheIndexList) :-
	getExampleIDs(TheAAList,TheIndexList),
	ppAAList(TheAAList).


